classdef ContinuousModel < handle
    % ContinuousModel  Wraps xdot = f(x,u) with an integrator.
    %
    %   model = ContinuousModel(f, nx, nu, integrator)
    %
    %   The integrator must expose step(f, x, u) -> x_next (row vectors).
    %   RK4 is the default one shipped.

    properties (SetAccess = private)
        f
        nx
        nu
        integrator
    end

    methods
        function obj = ContinuousModel(f, nx, nu, integrator)
            arguments
                f          function_handle
                nx   (1,1) double {mustBeInteger, mustBePositive}
                nu   (1,1) double {mustBeInteger, mustBePositive}
                integrator
            end
            assert(ismethod(integrator, 'step'), ...
                'ContinuousModel:integrator', ...
                'Integrator must implement step(f, x, u).');
            obj.f          = f;
            obj.nx         = nx;
            obj.nu         = nu;
            obj.integrator = integrator;
        end

        function x_next = step(obj, x, u)
            x_next = obj.integrator.step(obj.f, x, u);
        end
    end
end
