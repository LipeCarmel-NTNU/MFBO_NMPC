classdef InputMovementCost < handle
    % InputMovementCost  Σ‖Δu(k)‖_S² with u(0) = u_prev.

    properties (SetAccess = private)
        S
    end

    methods
        function obj = InputMovementCost(S)
            arguments
                S double
            end
            obj.S = S;
        end

        function J = evaluate(obj, parts, ctx)
            du = diff([ctx.u_prev; parts.u], [], 1);
            J  = 0;
            for k = 1 : size(du, 1)
                d = du(k, :).';
                J = J + d.' * obj.S * d;
            end
        end
    end
end
