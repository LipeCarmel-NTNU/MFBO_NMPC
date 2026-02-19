classdef NMPC_refactor_terminal < NMPC_refactor_abstract
    properties
        %% Set-points
        xsp = [1 5 1 0]
        usp = zeros(1, 3)

        %% Tuning parameters
        p = 60
        m = 6

        Q   = diag([10 1 1])
        Ru  = diag([2 2 1])
        Rdu = diag([100, 100 10])

        %% System parameters
        Ts = 1/60
        model

        nx = 4
        nu = 3

        %% Constraints
        Ymin = [0.5 0 -0.1]
        Ymax = [2 50 20]

        umax = 0.4 * [1 1 1]
        umin = zeros(1, 3)

        %% Optimizer
        optimizer_options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', ...
                'MaxFunEvals', Inf, 'MaxIterations', 1000, ...
                'StepTolerance', 1e-9, ...
                'OptimalityTolerance', 1e-6, ...
                'ScaleProblem', true)

        %% Debug
        debug_x = [1 8 1]
        debug_u = [1e-3 1e-3 2.1e-3]

        %% Terminal cost matrix
        P

        %% Scaling
        x_scale = [1 20 1]
        u_scale = [0.4 0.4 0.4]
    end

    methods
        function obj = NMPC_refactor_terminal(model, nx, nu)
            obj = obj.init(model, nx, nu);
        end

        function L = objfun(obj, w, u_init)
            x = reshape(w(1:obj.nx * (obj.p + 1)), [], obj.nx);
            u = reshape(w(obj.nx * (obj.p + 1) + 1:end), [], obj.nu);
            delta_u = diff([u_init; u], [], 1);

            L = 0;

            for k = 0:obj.p-1
                e = x(k+1,:) - obj.xsp;
                L = L + e * obj.Q * e.';
            end

            z = [x(end, :) - obj.xsp, u(end, :) - obj.usp];
            L = L + z * obj.P * z.';

            for k = 1:obj.m
                uk  = (u(k,:) - obj.usp).';
                duk = delta_u(k,:).';
                L = L + uk.' * obj.Ru * uk + duk.' * obj.Rdu * duk;
            end
        end
    end
end