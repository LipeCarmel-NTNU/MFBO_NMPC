classdef RK4 < handle
    % RK4  Fixed-step Runge-Kutta 4 integrator.

    properties (SetAccess = private)
        Ts
    end

    methods
        function obj = RK4(Ts)
            arguments
                Ts (1,1) double {mustBePositive, mustBeFinite}
            end
            obj.Ts = Ts;
        end

        function x_next = step(obj, f, x, u)
            % f : @(x,u) -> xdot, where x is a 1×nx row
            h = obj.Ts;
            k1 = RK4.row(f, x,            u);
            k2 = RK4.row(f, x + 0.5*h*k1, u);
            k3 = RK4.row(f, x + 0.5*h*k2, u);
            k4 = RK4.row(f, x +     h*k3, u);
            x_next = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end

    methods (Static, Access = private)
        function v = row(f, x, u)
            v = reshape(f(x, u), 1, []);
        end
    end
end
