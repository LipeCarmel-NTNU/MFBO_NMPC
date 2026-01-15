classdef MPC_abstract < handle
    properties (Abstract)
        %% Set-points
        % Depends on the formulation

        %% Tuning parameters
        p                         % Prediction horizon in steps
        m                         % Control  horizon in steps

        % Costs depend on the formulation

        %% System parameters
        Ts                         % Sampling time (h)
        model                      % System model in the form of a time invariant non-linear system: xdot = f(x, u)

        nx                         % Number of states
        nu                         % Number of manipulated variables
        
        %% Constraints
        Vmin
        Vmax

        Xmin
        Xmax

        % Numerical error may yield negligible negative sugar values. 
        % If this error becomes significant, reduce min_integrator_step, 
        % otherwise ignore it
        Smin
        Smax

        CO2min
        CO2max

        umax

        Lcon    % Lower bounds for the states (1 x nx)
        Ucon    % Upper bounds for the states (1 x nx)
        wL      % Lower bounds for the decision variables
        wU      % Upper bounds for the decision variables
        
        %% Integrator and Optimizer
        min_integrator_step    % Minimum integration step size required to solve the ode (h)

        optimizer_options

        latest_wopt

        %% DEBUG
        debug_x     % States to use for debugging
        debug_u     % Inputs to use for debugging
    end

    properties
        latest_flag = NaN
    end

    properties (Dependent)
        tf                              % Final time (h)
        n_integrator_steps              % Number of integration steps per multiple shooting subinterval
    end

    methods (Abstract)
        objfun(obj)
    end
    
    methods

        function [uk, x, u] = solve(obj, x_init, u_init)
            %
            % [uk, x, u] = obj.solve(x_init, u_init)
            % This method computes the optimal control inputs (uk) and the corresponding state trajectory (x) 
            % and control inputs (u) based on the initial state (x_init) and initial control inputs (u_init).
            %
            % Inputs:
            %   x_init - Initial state of the system (1 x nx vector), i.e: x(k|k)
            %   u_init - Initial control inputs (1 x nu vector) i.e: u(k-1)
            %
            % Outputs:
            %   uk    - Optimal control inputs for the next time step
            %   x     - State trajectory over the prediction horizon
            %   u     - Control inputs over the prediction horizon
            %
            % The method utilizes a moving horizon approach, where the 
            % previous solution is used as a guess 
            % for the optimization problem.
            %
            % The method also handles cases of infeasibility by first
            % satisfying continuity constraints then trying to solve the
            % optimization problem again. If the optimization fails again a
            % warning is thrown and u(k+1|k-1) is returned. Consider
            % adjusting the optimizer_options if such cases (mainly 
            % increasing MaxIter) 

            
            % Use the previous solution as a guess if it is available
            if isempty(obj.latest_wopt)
                % Solve the first iteration with no information
                [uk, x, u] = solve_first_iter(obj, x_init, u_init);

            elseif obj.latest_flag <= 0
                % If there is no previous solution, solve withuot
                % information
                failed_before = true;
                [uk, x, u] = solve_first_iter(obj, x_init, u_init);

            else
                % If there is a feasible past solution, apply moving
                % horizon

                % Get the previous u
                [x_pred, u] = data_from_w(obj, obj.latest_wopt);

                % Moving horizon
                u = [u(2:end, :); u(end, :)]; % if m=1 u(2:end, :) return a 0 x nu empty vector

                % Ensure continuity
                w0 = guess_from_initial(obj, x_init, u);
                [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init);
            end
            
            % In case of infeasibility, solve the ode with the u from the 
            % current solution thus ensuring satisfaction of the continuity
            % constraints. This tends to guide the optimizer to a feasible
            % solution
            if obj.latest_flag == -2 % < 1
                % Ensure continuity
                w0 = guess_from_initial(obj, x_init, u);
                % Retry
                [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init);
            end
            if obj.latest_flag < 1
                if failed_before
                    warning('MPC - FAILED TWICE')
                    uk = zeros(1, obj.nu);
                else
                    warning('MPC - Failed to solve. Using the u(k+1|k-1)')
                    % u = u_init;
                    [~, u] = data_from_w(obj, obj.latest_wopt);
                    if obj.m > 1
                        uk = u(2, :);
                    else
                        uk = u(1, :);
                    end
                end
            end

        end
        
        function [uk, x, u, fval] = solve_first_iter(obj, x_init, u_init)
            % Solve
            obj.constraints();
            w0 = guess_from_initial(obj, x_init, u_init);
            [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init);

        end
        
        function [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init)
            % is_infeasible(obj, w0, x_init)
            [wopt,fval,exitflag] = fmincon(@(w) obj.objfun(w, u_init),w0,...
                    [],[],[],[],obj.wL,obj.wU,...
                    @(w) obj.confun(w, x_init),...
                    obj.optimizer_options);
            if exitflag >= 0
                obj.latest_wopt = wopt;
            end
            obj.latest_flag = exitflag;
            
            % Extract states and controls
            [x, u] = obj.data_from_w(wopt);
            obj.validate(false, x, u);
            uk = u(1, :);
        end
        
        function tf = get.tf(obj)
            % GET.TF Getter for final time
            %   Computes the final time based on the prediction horizon and
            % the sampling time.
            %
            %   Outputs:
            %       tf - Final time  (h)

            tf = obj.p * obj.Ts;
        end
        
        function n = get.n_integrator_steps(obj)
            % Computes the number of integrator steps to preserve numerical
            % precision.
            n = max(1, round(obj.Ts / obj.min_integrator_step));
        end
        
        function [c, ceq] = confun(obj, w, x_init, u_init)
            % CONFUN Constraint function
            %   Computes the continuity and initial state constraints for 
            % the optimization problem (compatible with fmincon).
            %
            %   Inputs:
            %       w - Decision variable vector
            %       x_init - Initial state vector
            %
            %   Outputs:
            %       c   -  Non-linear inequality constraints (empty)
            %       ceq - Non-linear equality constraints

            % Extract states and control actions
            [x, u] = obj.data_from_w(w);
        
            
            % Initialize states at nodes
            states_atNodes = zeros(obj.p + 1, obj.nx);
            
            % Single shooting in each subinterval using RK4 integration
            for i = 1:obj.p
                x0 = x(i, :);
                if i <= obj.m
                    uk = u(i, :);
                end
                
                % Integration between the nodes
                states = obj.integrator(x0, uk);

                states_atNodes(i + 1, :) = states(end, :);
            end
            
            % Continuity constraints
            ceq_temp = x(2:end, :) - states_atNodes(2:end, :);
            ceq_temp = [ceq_temp; x(1, :) - x_init];
            
            % Equality constraints
            ceq = reshape(ceq_temp, [], 1);
            
            % Non-linear inequality constraints (none)
            c = [];
        end

        
        function [states_atNodes] = confun_initial(obj, x_init, u)
            % CONFUN_INITIAL Initial constraint function
            %   Computes the states at nodes based on initial conditions for the optimization problem.
            %
            %   Inputs:
            %       x_init - Initial state vector
            %       u - Control matrix (m x nu)
            %
            %   Outputs:
            %       states_atNodes - States at each node
            
            % Initialize states at nodes
            states_atNodes = [x_init; zeros(obj.p, obj.nx)];
            
            % Single shooting in each subinterval using RK4 integration
            for i = 1:obj.p
                x0 = states_atNodes(i, :);
                if i <= obj.m
                    uk = u(i, :);
                end
                
                % Integration between the nodes
                states = obj.integrator(x0, uk);

                states_atNodes(i + 1, :) = states(end, :);
            end
        end

        function [x, u] = data_from_w(obj, w)
            % Returns [x, u] from the decision variable vector.

            len_x =  obj.nx * (obj.p + 1);  % length of x vector (nx initial states + p predictions)
            len_u = obj.m*obj.nu;           % length of u vector (nu * m control actions)

            x = reshape(w(1 : len_x), [], obj.nx);
            u = reshape(w(len_x + 1: len_x + len_u), [], obj.nu);
        end

        function [wL, wU] = constraints(obj)
            % CONSTRAINTS Set lower bound and upper bound for states and control
            %   Computes the lower and upper contraints for the decision variables.
            %
            %   Outputs:
            %       wL - Lower contraints vector for the decision variables
            %       wU - Upper contraints vector for the decision variables

            obj.latest_wopt = []; % incompatible dimensions

            % Dimensions of decision variables
            len_x =  obj.nx * (obj.p + 1);  % length of x vector (nx initial states + p predictions)
            len_u = obj.m*obj.nu;           % length of u vector (nu * m control actions)

            %% Lower and upper bouds
            
            % Lower bound for states
            xL =  ones(obj.p + 1, obj.nx).*obj.Lcon;
            xL = reshape(xL ,[], 1);
            
            % Upper bound for states
            xU =  ones(obj.p + 1, obj.nx).*obj.Ucon;
            xU = reshape(xU ,[], 1);
            
            
            % Lower bound for control actions
            uL =   0 * ones(len_u, 1);

            % Upper bound for control actions
            uU = ones(obj.m, obj.nu).*obj.umax;
            uU = reshape(uU ,[], 1);      

            wL = [xL; uL];
            wU = [xU; uU];
            
            obj.wL = wL;
            obj.wU = wU;

            %% Linear inequality constraints
            % Aeq = [zeros(len_u, len_x), eye(len_u), - eye(len_u)]; % 0 + u - s <= umax
            % beq = reshape(ones(obj.m, obj.nu).*obj.umax, [], 1);
            % obj.Aeq = Aeq;
            % obj.beq = beq;


        end

        function w0 = guess_from_initial(obj, x_init, u_init)
            % GUESS Initial guess for the decision variables
            %   Provides an initial guess for the optimization problem based on initial state and control.
            %
            %   Inputs:
            %       x_init - Initial state vector
            %       u_init - Control action vector at k-1, or m x nu matrix
            %
            %   Outputs:
            %       w0 - Initial guess vector for the decision variables

            if size(u_init, 1) == 1
                u = [u_init; zeros(obj.m - 1, obj.nu)];                  % Initial guess for control actions
            else
                u = u_init;
            end
            u(u < 0) = 0;
            [states_atNodes] = obj.confun_initial(x_init, u);        % Compute states at nodes
            w0_x = reshape(states_atNodes, [], 1);                   % Reshape states to a vector
            w0_u = reshape(u, [], 1);                                % Reshape control actions to a vector
            w0 = [w0_x; w0_u];                                       % Combine states and control actions for initial guess
        end
        
        function L = eval_cost(obj, x_init, u_init, u)
            arguments
                obj 
                x_init 
                u_init 
                u = []
            end
            if isempty(u)
                w = guess_from_initial(obj, x_init, u_init);
            else
                w = guess_from_initial(obj, x_init, u);
            end
            [L] = objfun(obj, w, u_init);
        end
        function test_dimensions(obj)
            assert(length(obj.debug_u) == obj.nu, sprintf('Incorrect dimensions for debug_u (%d) and/or nu (%d)', length(obj.debug_u), obj.nu));
            assert(length(obj.debug_x) == obj.nx, sprintf('Incorrect dimensions for debug_x (%d) and/or nx (%d)', length(obj.debug_x), obj.nx));
            assert(length(obj.Lcon) == obj.nx, sprintf('Incorrect dimensions for Lcon (%d) and/or nx (%d)', length(obj.Lcon), obj.nx));
            assert(length(obj.Ucon) == obj.nx, sprintf('Incorrect dimensions for Ucon (%d) and/or nx (%d)', length(obj.Ucon), obj.nx));
            assert(length(obj.umax) == obj.nu, sprintf('Incorrect dimensions for umax (%d) and/or nu (%d)', length(obj.umax), obj.nu));

        end
        
        function test_guess(obj)
            % Solve the open loop problem
            w0 = guess_from_initial(obj, obj.debug_x, obj.debug_u);

            % Calculate continuity constraints
            [c, ceq] = confun(obj, w0, obj.debug_x);

            assert(all(ceq == 0), "Non-linear constraints not satisfied for the guess solution")

        end

        function success = validate(obj, quiet, x, u)
            success = true;             % Flag indicating successful integration

            if nargin == 1
                quiet = false;          % Bool for printing
            end

            % Solve the open loop problem
            if nargin < 3
                % x and u were not provided, run the default test case
                w0 = guess_from_initial(obj, obj.debug_x, obj.debug_u);
            else
                w0 = guess_from_initial(obj, x(1, :), u);
            end

            % Extract data
            [x, u] = data_from_w(obj, w0);

            % Solve again with a built-in solver
            x_ode = x;
            for i = 1 : obj.p
                % Current state is the initial condition
                x0 = x_ode(i, :);
                if i <= obj.m
                    uk = u(i,:);
                end
                [~, y] = ode45(@(t,x) obj.model(x,uk), [0 obj.Ts], x0);
                % Save next state (next initial condition)
                x_ode(i+1, :) = y(end,:);
            end
            e = x - x_ode;

            % Indexes
            max_e = max(abs(e(end,:)));
            mae = mean(abs(e));

            if not(quiet)
                fprintf('\n\n')
                disp('Validating MPC predictions')
                disp('')
                disp('Maximum Absolute Error:')
                disp(max_e)
                disp('Mean Absolute Error:')
                disp(mae)
            end

            % Warn in case of error > tol
            tol = 0.001; % A fraction of Ks
            if any(max_e > tol)
                success = false;
                warning('DEBUGGING IS NECESSARY')
            end

        end
   
        function states = integrator(obj, x0, uk)
            h = obj.Ts;
            
            states = zeros(N + 1, obj.nx);
            states(1, :) = x0;

            [~,y] = RKF45_book(@(t,x) obj.model(x, uk), [0 h], x0, h);

            states(end,:) = y(end,:);
        end

        function [infeasible] = is_infeasible(obj, w0, x_init)
            function check(violation_bool, zero_vec)
                error_idx = find(violation_bool);
                worst = max(abs(zero_vec(error_idx)));
                if worst > 0
                    warning('High constraint violation')
                    disp(worst)
                end
            end
            infeasible = false;
            [~, ceq] = obj.confun(w0, x_init);

            ceq_violation = ceq ~= 0;
            check(ceq_violation, ceq)

            wL_violation = w0 < obj.wL;
            check(wL_violation, w0 - obj.wL)

            wU_violation = obj.wU < w0;
            check(wU_violation, w0 - obj.wU)


            infeasible = bitor(infeasible, any(ceq_violation));
            infeasible = bitor(infeasible, any(wL_violation));
            infeasible = bitor(infeasible, any(wU_violation));
            if infeasible
                keyboard
            end
        end


    end
    methods (Static)
        function J = norm_Q(x, Q)
            % NORM_Q Quadratic norm calculation
            %   Computes the quadratic norm of a vector with a given weight.
            %
            %   Inputs:
            %       x - Input vector
            %       Q - Weight matrix/scalar
            %
            %   Outputs:
            %       J - Quadratic norm value

            J = x' * Q * x;
        end
        
        function w = w_from_data(x, u)
            w_x = reshape(x, [], 1);                   % Reshape states to a vector
            w_u = reshape(u, [], 1);                   % Reshape control actions to a vector
            w = [w_x; w_u];  
        end
    end
end

