classdef NMPC_mono_dilution < NMPC_abstract
    properties
        %% Set-points
        xsp = [1 5 1 0]                % Setpoint for states (1 x nx)

        %% Tuning parameters
        p = 60;                         % Prediction horizon in steps
        m = 6;                          % Control  horizon in steps

        Q   = diag([10 1 1 0])          % State tracking weights (nx x nx)
        Rdu = diag([10, 0.1 0.1]);      % Input increment weights (nu x nu)

        %% System parameters
        Ts = 1/60;                      % Sampling time (h)
        model                           % System model in the form of a time invariant non-linear system: xdot = f(x, u)

        %% Constraints

        % Numerical error may yield negligible negative sugar values.
        % If this error becomes significant, reduce min_integrator_step,
        % otherwise ignore it

        Ymin = [0.5 0 -0.1 0]
        Ymax = [2 50 20 Inf]

        umax = 0.4*[1 1 1]
        umin = zeros(1, 3)
        
        dilution = 1e7

        wL      % Lower bounds for the decision variables
        wU      % Upper bounds for the decision variables

        %% Integrator and Optimizer
        min_integrator_step = 0.007;    % Minimum integration step size required to solve the ode (h)

        optimizer_options = optimoptions('fmincon','Display','Iter','Algorithm','sqp',...
                'MaxFunEvals',Inf, 'MaxIterations', 100);

        latest_wopt

        %% DEBUG
        debug_x = [1 8 1 0.04]          % States to use for debugging
        debug_u = [1e-3 1e-3 2.1e-3]    % Inputs to use for debugging

    end
    methods
        function obj = NMPC_mono_dilution(model, nx, nu)
            % Initialization inherited from abstract class
            obj = obj.init(model, nx, nu);
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

            % Sum form
            L = 0;
            for k = 1 : obj.p + 1
                e = x(k,:) - obj.xsp;
                L = L + e * obj.Q * e.';
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
