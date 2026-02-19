classdef NMPC_refactor_abstract < handle
    properties (Abstract)
        % Tuning parameters
        p
        m

        % System parameters
        Ts
        model
        nx
        nu

        % Constraints
        Ymin
        Ymax
        umin
        umax

        % Optimizer
        optimizer_options

        % Debug values
        debug_x
        debug_u
    end

    properties
        % Tracking references and weights (set by concrete class)
        xsp
        usp
        Q
        Ru
        Rdu

        % Scaling (physical -> scaled: v_scaled = v ./ v_scale)
        x_scale = []
        u_scale = []
        scale_constraints = true

        % State
        latest_wopt = []
        latest_flag = NaN
        real_time_validation = false
        quiet_validation = false

        % Bounds in scaled coordinates
        wL
        wU
    end

    properties (Dependent)
        tf
    end

    methods (Abstract)
        L = objfun(obj, w, u_init)
    end

    methods
        function obj = init(obj, model, nx, nu)
            obj.model = model;
            obj.nx = nx;
            obj.nu = nu;

            if isempty(obj.x_scale)
                obj.x_scale = ones(1, obj.nx);
            end
            if isempty(obj.u_scale)
                obj.u_scale = ones(1, obj.nu);
            end
            obj.validate_scales();

            obj.rebuild_bounds();
            obj.test_dimensions();
            obj.validate();
        end

        function [uk, x, u] = solve(obj, x_init, u_init)
            failed_before = (obj.latest_flag < 0);

            if isempty(obj.latest_wopt) || failed_before
                [uk, x, u] = obj.solve_first_iter(x_init, u_init);
            else
                [~, u] = obj.data_from_w(obj.latest_wopt);
                u = [u(2:end, :); u(end, :)];
                w0 = obj.guess_from_initial(x_init, u);
                [uk, x, u] = obj.solve_optimization(w0, x_init, u_init);
            end

            if obj.latest_flag == -2
                w0 = obj.guess_from_initial(x_init, u);
                [uk, x, u] = obj.solve_optimization(w0, x_init, u_init);
            end

            if obj.latest_flag < 0
                if failed_before
                    warning('MPC - FAILED AGAIN')
                    uk = zeros(1, obj.nu);
                else
                    warning('MPC - Failed to solve. Using the u(k+1|k-1)')
                    if isempty(obj.latest_wopt)
                        uk = zeros(1, obj.nu);
                    else
                        [~, u] = obj.data_from_w(obj.latest_wopt);
                        if obj.m > 1
                            uk = u(2, :);
                        else
                            uk = u(1, :);
                        end
                    end
                end
            end
        end

        function [uk, x, u, fval] = solve_first_iter(obj, x_init, u_init)
            obj.rebuild_bounds();
            w0 = obj.guess_from_initial(x_init, u_init);
            [uk, x, u, fval] = obj.solve_optimization(w0, x_init, u_init);
        end

        function [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init)
            w0s = obj.scale_w(w0);

            [wopts, fval, exitflag] = fmincon(
                @(ws) obj.objfun_scaled(ws, u_init), w0s, ...
                [],[],[],[], obj.wL, obj.wU, ...
                @(ws) obj.confun_scaled(ws, x_init), ...
                obj.optimizer_options);

            wopt = obj.unscale_w(wopts);
            if exitflag >= 0
                obj.latest_wopt = wopt;
            end
            obj.latest_flag = exitflag;

            [x, u] = obj.data_from_w(wopt);
            if obj.real_time_validation
                obj.validate(obj.quiet_validation, x, u);
            end
            uk = u(1, :);
        end

        function [c, ceq] = confun(obj, w, x_init)
            [x, u] = obj.data_from_w(w);
            states_atNodes = zeros(obj.p + 1, obj.nx);

            for i = 1:obj.p
                x0 = x(i, :);
                if i <= obj.m
                    uk = u(i, :);
                end
                states = obj.integrator(x0, uk);
                states_atNodes(i + 1, :) = states(end, :);
            end

            ceq_temp = x(2:end, :) - states_atNodes(2:end, :);
            ceq_temp = [ceq_temp; x(1, :) - x_init];

            ceq = reshape(ceq_temp, [], 1);
            c = [];
        end

        function [c, ceq] = confun_scaled(obj, ws, x_init)
            w = obj.unscale_w(ws);
            [c_phys, ceq_phys] = obj.confun(w, x_init);

            if obj.scale_constraints
                c = c_phys;
                ceq = obj.scale_ceq(ceq_phys);
            else
                c = c_phys;
                ceq = ceq_phys;
            end
        end

        function J = objfun_scaled(obj, ws, u_init)
            w = obj.unscale_w(ws);
            J = obj.objfun(w, u_init);
        end

        function [wL, wU] = constraints(obj)
            len_x = obj.nx * (obj.p + 1);
            len_u = obj.m * obj.nu;

            xL = ones(obj.p + 1, obj.nx) .* obj.Ymin;
            xU = ones(obj.p + 1, obj.nx) .* obj.Ymax;
            uL = ones(obj.m, obj.nu) .* obj.umin;
            uU = ones(obj.m, obj.nu) .* obj.umax;

            wL = [reshape(xL, [], 1); reshape(uL, [], 1)];
            wU = [reshape(xU, [], 1); reshape(uU, [], 1)];

            obj.wL = obj.scale_w(wL);
            obj.wU = obj.scale_w(wU);
        end

        function rebuild_bounds(obj)
            obj.latest_wopt = [];
            obj.constraints();
        end

        function [states_atNodes] = confun_initial(obj, x_init, u)
            states_atNodes = [x_init; zeros(obj.p, obj.nx)];
            for i = 1:obj.p
                x0 = states_atNodes(i, :);
                if i <= obj.m
                    uk = u(i, :);
                end
                states = obj.integrator(x0, uk);
                states_atNodes(i + 1, :) = states(end, :);
            end
        end

        function x_next = integrator(obj, x0, uk)
            h = obj.Ts;
            [~, y] = RKF45_book(@(t, x) obj.model(x, uk), [0 h], x0, h);
            x_next = y(end, :);
        end

        function [x, u] = data_from_w(obj, w)
            len_x = obj.nx * (obj.p + 1);
            len_u = obj.m * obj.nu;
            x = reshape(w(1:len_x), [], obj.nx);
            u = reshape(w(len_x + 1:len_x + len_u), [], obj.nu);
        end

        function w = w_from_data(obj, x, u)
            w = [reshape(x, [], 1); reshape(u, [], 1)];
        end

        function w0 = guess_from_initial(obj, x_init, u_init)
            if size(u_init, 1) == 1
                u = [u_init; zeros(obj.m - 1, obj.nu)];
            else
                u = u_init;
            end
            u(u < 0) = 0;
            states_atNodes = obj.confun_initial(x_init, u);
            w0_x = reshape(states_atNodes, [], 1);
            w0_u = reshape(u, [], 1);
            w0 = [w0_x; w0_u];
        end

        function success = validate(obj, quiet, x, u)
            success = true;
            if nargin == 1
                quiet = false;
            end

            if nargin < 3
                w0 = obj.guess_from_initial(obj.debug_x, obj.debug_u);
            else
                w0 = obj.guess_from_initial(x(1, :), u);
            end

            [x, u] = obj.data_from_w(w0);
            opts = odeset('Nonnegative', 1:obj.nx);
            x_ode = x;
            x_ode(x_ode < 0) = 0;
            for i = 1:obj.p
                x0 = x_ode(i, :);
                if i <= obj.m
                    uk = u(i, :);
                end
                [~, y] = ode45(@(t, x) obj.model(x, uk), [0 obj.Ts], x0, opts);
                x_ode(i + 1, :) = y(end, :);
            end
            e = x - x_ode;
            max_e = max(abs(e(end, :)));
            mae = mean(abs(e));

            if ~quiet
                fprintf('\n\n')
                disp('Validating MPC predictions')
                disp('')
                disp('Maximum Absolute Error:')
                disp(max_e)
                disp('Mean Absolute Error:')
                disp(mae)
            end

            tol = 0.001;
            if any(max_e > tol)
                success = false;
                warning('DEBUGGING IS NECESSARY')
            end
        end

        function test_dimensions(obj)
            assert(length(obj.debug_u) == obj.nu, sprintf('Incorrect dimensions for debug_u (%d) and/or nu (%d)', length(obj.debug_u), obj.nu));
            assert(length(obj.debug_x) == obj.nx, sprintf('Incorrect dimensions for debug_x (%d) and/or nx (%d)', length(obj.debug_x), obj.nx));
            assert(length(obj.Ymin) == obj.nx, sprintf('Incorrect dimensions for Ymin (%d) and/or nx (%d)', length(obj.Ymin), obj.nx));
            assert(length(obj.Ymax) == obj.nx, sprintf('Incorrect dimensions for Ymax (%d) and/or nx (%d)', length(obj.Ymax), obj.nx));
            assert(length(obj.umax) == obj.nu, sprintf('Incorrect dimensions for umax (%d) and/or nu (%d)', length(obj.umax), obj.nu));
        end

        function validate_scales(obj)
            assert(isvector(obj.x_scale) && numel(obj.x_scale) == obj.nx, ...
                'x_scale must be (1 x nx).');
            assert(isvector(obj.u_scale) && numel(obj.u_scale) == obj.nu, ...
                'u_scale must be (1 x nu).');
            assert(all(isfinite(obj.x_scale)) && all(obj.x_scale > 0), ...
                'x_scale must be > 0 and finite.');
            assert(all(isfinite(obj.u_scale)) && all(obj.u_scale > 0), ...
                'u_scale must be > 0 and finite.');
        end

        function xs = scale_x(obj, x)
            xs = x ./ obj.x_scale;
        end

        function x = unscale_x(obj, xs)
            x = xs .* obj.x_scale;
        end

        function us = scale_u(obj, u)
            us = u ./ obj.u_scale;
        end

        function u = unscale_u(obj, us)
            u = us .* obj.u_scale;
        end

        function ws = scale_w(obj, w)
            [x, u] = obj.data_from_w(w);
            xs = obj.scale_x(x);
            us = obj.scale_u(u);
            ws = obj.w_from_data(xs, us);
        end

        function w = unscale_w(obj, ws)
            [xs, us] = obj.data_from_w(ws);
            x = obj.unscale_x(xs);
            u = obj.unscale_u(us);
            w = obj.w_from_data(x, u);
        end

        function ceqs = scale_ceq(obj, ceq)
            if isempty(ceq)
                ceqs = ceq;
                return
            end\n            s = obj.ceq_scale_vector();\n            ceqs = ceq ./ s;\n        end\n\n        function s = ceq_scale_vector(obj)\n            sx = obj.x_scale(:);\n            s_cont = repmat(sx, obj.p, 1);\n            s_init = sx;\n            s = [s_cont; s_init];\n        end\n\n        function tf = get.tf(obj)\n            tf = obj.p * obj.Ts;\n        end\n    end\n\n    methods (Static)\n        function J = norm_Q(x, Q)\n            J = x' * Q * x;\n        end\n    end\nend
