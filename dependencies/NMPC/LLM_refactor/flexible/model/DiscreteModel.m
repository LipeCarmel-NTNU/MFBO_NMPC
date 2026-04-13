classdef DiscreteModel < handle
    % DiscreteModel  x_next = f(x,u); no integrator.

    properties (SetAccess = private)
        f
        nx
        nu
    end

    methods
        function obj = DiscreteModel(f, nx, nu)
            arguments
                f          function_handle
                nx   (1,1) double {mustBeInteger, mustBePositive}
                nu   (1,1) double {mustBeInteger, mustBePositive}
            end
            obj.f  = f;
            obj.nx = nx;
            obj.nu = nu;
        end

        function x_next = step(obj, x, u)
            x_next = reshape(obj.f(x, u), 1, []);
        end
    end
end
