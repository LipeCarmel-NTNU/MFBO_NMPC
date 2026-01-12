classdef MPC_mono_dilution < MPC_abstract
    properties
        %% Set-points
        Vsp = 1                        % Setpoint for volume (L)
        Xsp = 5                        % Setpoint for biomass (g/L)
        Ssp = 1                        % Setpoint for biomass (g/L)

        %% Tuning parameters
        p = 60;                         % Prediction horizon in steps
        m = 6;                          % Control  horizon in steps

        Q_V = 10
        Q_X = 1
        Q_S = 1
        Rdu = diag([10, 0.1 0.1]);

        %% System parameters
        Ts = 1/60;                     % Sampling time (h)
        model                          % System model in the form of a time invariant non-linear system: xdot = f(x, u)

        nx = 4;                        % Number of states
        nu = 3;                        % Number of manipulated variables
        
        %% Constraints
        Vmin = 0.5
        Vmax = 2

        Xmin = 0
        Xmax = 50

        % Numerical error may yield negligible negative sugar values. 
        % If this error becomes significant, reduce min_integrator_step, 
        % otherwise ignore it
        Smin = -0.1
        Smax = 20

        CO2min = 0
        CO2max = Inf

        umax = 0.4*[1 1 1]

        dilution = 1e7

        Lcon    % Lower bounds for the states (1 x nx)
        Ucon    % Upper bounds for the state (1 x nx)
        wL      % Lower bounds for the decision variables
        wU      % Upper bounds for the decision variables
        
        %% Integrator and Optimizer
        min_integrator_step = 0.007;    % Minimum integration step size required to solve the ode (h)

        optimizer_options = optimoptions('fmincon','Display','Iter','Algorithm','sqp',...
                'MaxFunEvals',Inf, 'MaxIterations', 100);

        latest_wopt

        %% DEBUG
        debug_x = [1 8 1 0.04]      % States to use for debugging
        debug_u = [1e-3 1e-3 2.1e-3]     % Inputs to use for debugging

    end
    methods
        function obj = MPC_mono_dilution(model)
            % MPC Constructor
            %   Initializes the MPC object with the provided model and position setpoints.
            %
            %   Inputs:
            %       model - System model function handle
            %       Vsp - Setpoint for volume (L)
            %       Xsp - Setpoint for biomass (g/L)x

            obj.model = model;

            % Upper and lower constraints
            obj.Lcon = [obj.Vmin, obj.Xmin, obj.Smin, obj.CO2min];
            obj.Ucon = [obj.Vmax, obj.Xmax, obj.Smax, obj.CO2max];

            [obj.wL, obj.wU] = obj.constraints();


            obj.test_dimensions()
            obj.validate();
        end
        function [L] = objfun(obj, w, u_init)
            % OBJFUN Objective function
            %   Computes the objective function value for the optimization problem.
            %
            %   Inputs:
            %       w - Decision variable vector
            %
            %   Outputs:
            %       L - Objective function value

            % Extract states and control actions
            x = reshape(w(1:obj.nx * (obj.p + 1)), [], obj.nx);
            u = reshape(w(obj.nx * (obj.p + 1) + 1:end), [], obj.nu);
            delta_u = diff([u_init; u],[],1);
            
            % Weights
            Qblk = blkdiag(obj.Q_V, obj.Q_X, obj.Q_S, 0);
            xsp  = [obj.Vsp, obj.Xsp, obj.Ssp, 0];
        
            % Sum form
            L = 0;
            for k = 1 : obj.p + 1
                e = x(k,:) - xsp;
                L = L + e * Qblk * e.';
            end
            for k = 1:obj.m
                du = delta_u(k,:).';
                L = L + du.' * obj.Rdu * du;
            end
        
            % Slacked dilution
            L = L + obj.dilution * sum(u(:,2));
        end
    end
end