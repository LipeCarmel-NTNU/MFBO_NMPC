classdef EconomicCost < handle
    % EconomicCost  Σ w · φ(x_k, u_k) over the prediction horizon.
    %
    %   c = EconomicCost(phi)                     % scalar-valued @(x,u) -> R
    %   c = EconomicCost(phi, weight=w)           % scalar weight (default 1)

    properties (SetAccess = private)
        phi
        weight
    end

    methods
        function obj = EconomicCost(phi, args)
            arguments
                phi        function_handle
                args.weight (1,1) double = 1
            end
            obj.phi    = phi;
            obj.weight = args.weight;
        end

        function J = evaluate(obj, parts, ctx)
            x = parts.x;
            u = parts.u;
            p = ctx.p;
            mh = ctx.m;
            J  = 0;
            uk = u(1, :);
            for k = 1 : p
                if k <= mh
                    uk = u(k, :);
                end
                J = J + obj.weight * obj.phi(x(k, :), uk);
            end
        end
    end
end
