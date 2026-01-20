classdef NMPC_terminal < NMPC_abstract_scaled
    properties
        %% Set-points
        xsp = [1 5 1 0]                 % Setpoint for states (1 x nx)
        usp = zeros(1,3)                % Reference for input penalty (1 x nu)

        %% Tuning parameters
        p = 60;                         % Prediction horizon in steps
        m = 6;                          % Control  horizon in steps

        Q   = diag([10 1 1])             % State tracking weights (nx x nx)
        Ru  = diag([2 2 1]);             % Input weights (nu x nu)
        Rdu = diag([100, 100 10]);       % Input increment weights (nu x nu)

        %% System parameters
        Ts = 1/60;                       % Sampling time (h)
        model

        nx = 4;
        nu = 3;

        %% Constraints
        Ymin = [0.5 0 -0.1]
        Ymax = [2 50 20]

        umax = 0.4*[1 1 1]
        umin = zeros(1, 3)

        wL
        wU

        %% Integrator and Optimizer

        optimizer_options = optimoptions('fmincon','Display','Iter','Algorithm','sqp', ...
                'MaxFunEvals',Inf, 'MaxIterations', 1000, ...
                'StepTolerance', 1e-9, ...
                'OptimalityTolerance',1e-6,...
                'ScaleProblem', true);

        latest_wopt

        %% DEBUG
        debug_x = [1 8 1]
        debug_u = [1e-3 1e-3 2.1e-3]

        %% Terminal cost matrix
        P
    end
    methods
        function obj = NMPC_terminal(model, nx, nu)
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
            delta_u = diff([u_init; u], [], 1);

            L = 0;

            % Tracking cost
            for k = 0 : obj.p - 1
                e = x(k+1,:) - obj.xsp;
                L = L + e * obj.Q * e.';
            end
            % Terminal cost
            z = [x(end, :) - obj.xsp, u(end, :) - obj.usp];
            L = L + z * obj.P * z.';
            % Regularization
            for k = 1:obj.m
                uk  = (u(k,:) - obj.usp).';
                duk = delta_u(k,:).';
                L = L + uk.'*obj.Ru*uk + duk.'*obj.Rdu*duk;
            end

        end
    end
end
