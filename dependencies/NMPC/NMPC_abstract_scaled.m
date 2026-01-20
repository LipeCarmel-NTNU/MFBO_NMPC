classdef NMPC_abstract_scaled < NMPC_abstract
    properties
        % Elementwise scaling factors (physical -> scaled: v_scaled = v ./ v_scale)
        % Keep strictly positive.
        x_scale = [1 20 1]   % (1 x nx)
        u_scale = [0.4 0.4 0.4]   % (1 x nu)

        % If true, scale nonlinear equality constraints (recommended for SQP)
        scale_constraints = true
    end

    methods
        function obj = init(obj, model, nx, nu)
            obj = init@NMPC_abstract(obj, model, nx, nu);

            if isempty(obj.x_scale)
                obj.x_scale = ones(1, obj.nx);
            end
            if isempty(obj.u_scale)
                obj.u_scale = ones(1, obj.nu);
            end
            obj.validate_scales();

            % Rebuild bounds in scaled coordinates
            [obj.wL, obj.wU] = obj.constraints();
        end

        function [uk, x, u] = solve(obj, x_init, u_init)
            [uk, x, u] = solve@NMPC_abstract(obj, x_init, u_init);
        end

        function [uk, x, u, fval] = solve_optimization(obj, w0, x_init, u_init)
            w0s = obj.scale_w(w0);

            [wopts, fval, exitflag] = fmincon( ...
                @(ws) obj.objfun_scaled(ws, u_init), ...
                w0s, [],[],[],[], obj.wL, obj.wU, ...
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
            % Keep the physical (unscaled) constraint available for debugging.
            [c, ceq] = confun@NMPC_abstract(obj, w, x_init);
        end

        function [wL, wU] = constraints(obj)
            [wL_phys, wU_phys] = constraints@NMPC_abstract(obj);
            wL = obj.scale_w(wL_phys);
            wU = obj.scale_w(wU_phys);

            obj.wL = wL;
            obj.wU = wU;
        end
    end

    methods (Access = protected)
        function J = objfun_scaled(obj, ws, u_init)
            w = obj.unscale_w(ws);
            J = obj.objfun(w, u_init);
        end

        function [c, ceq] = confun_scaled(obj, ws, x_init)
            w = obj.unscale_w(ws);

            [c_phys, ceq_phys] = obj.confun(w, x_init);

            if obj.scale_constraints
                c   = obj.scale_c(c_phys);
                ceq = obj.scale_ceq(ceq_phys);
            else
                c   = c_phys;
                ceq = ceq_phys;
            end
        end
    end

    methods
        function validate_scales(obj)
            assert(isvector(obj.x_scale) && numel(obj.x_scale) == obj.nx, ...
                'x_scale must be (1 x nx).');
            assert(isvector(obj.u_scale) && numel(obj.u_scale) == obj.nu, ...
                'u_scale must be (1 x nu).');

            assert(all(isfinite(obj.x_scale)) && all(obj.x_scale > 0), 'x_scale must be > 0 and finite.');
            assert(all(isfinite(obj.u_scale)) && all(obj.u_scale > 0), 'u_scale must be > 0 and finite.');
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
            % ceq packs:
            %   [x(2:end,:) - xhat(2:end,:);
            %    x(1,:) - x_init]
            % All equality residuals are in state units.
            if isempty(ceq)
                ceqs = ceq;
                return
            end

            s = obj.ceq_scale_vector();
            ceqs = ceq ./ s;
        end

        function cs = scale_c(obj, c)
            % No nonlinear inequality constraints in base class.
            cs = c;
        end
    end

    methods (Access = private)
        function s = ceq_scale_vector(obj)
            % Build scaling for ceq (vectorised) based on x_scale.
            % Number of state equalities: (p * nx) for continuity + (nx) for initial condition.
            nx = obj.nx;
            p  = obj.p;

            sx = obj.x_scale(:);
            s_cont = repmat(sx, p, 1);
            s_init = sx;

            s = [s_cont; s_init];
        end
    end
end
