classdef QuadraticTracking < handle
    % QuadraticTracking  Σ‖y - y_sp‖_Q² over p + Σ‖u - u_sp‖_R² over m.
    %
    %   c = QuadraticTracking(Q, R, y_sp, u_sp)
    %   c = QuadraticTracking(..., h_y=@(x)...)    % output map, default y=x

    properties (SetAccess = private)
        Q
        R
        y_sp
        u_sp
        h_y
    end

    methods
        function obj = QuadraticTracking(Q, R, y_sp, u_sp, args)
            arguments
                Q  double
                R  double
                y_sp double
                u_sp double
                args.h_y function_handle = @(x) x
            end
            obj.Q    = Q;
            obj.R    = R;
            obj.y_sp = y_sp;
            obj.u_sp = u_sp;
            obj.h_y  = args.h_y;
        end

        function J = evaluate(obj, parts, ctx)
            x = parts.x;
            u = parts.u;
            p = ctx.p;
            mh = ctx.m;
            ny = ctx.ny;

            ysp = QuadraticTracking.expand(obj.y_sp, p,  ny);
            usp = QuadraticTracking.expand(obj.u_sp, mh, ctx.nu);

            J = 0;
            for k = 1 : p
                yk = obj.h_y(x(k, :));
                e  = (yk - ysp(k, :)).';
                J  = J + e.' * obj.Q * e;
            end
            for k = 1 : mh
                ek = (u(k, :) - usp(k, :)).';
                J  = J + ek.' * obj.R * ek;
            end
        end
    end

    methods (Static)
        function M = expand(sp, n_rows, n_cols)
            if size(sp, 1) == 1
                M = repmat(sp, n_rows, 1);
            else
                assert(size(sp, 1) == n_rows && size(sp, 2) == n_cols, ...
                    'QuadraticTracking:sp', 'Setpoint size mismatch.');
                M = sp;
            end
        end
    end
end
