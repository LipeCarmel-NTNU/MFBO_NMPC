classdef Continuity < handle
    % Continuity  Multiple-shooting continuity equalities.
    %
    %   Auto-added by the orchestrator. Residuals are divided by x_scale
    %   so they are dimensionless (matches NMPC_abstract_scaled).

    methods
        function [c, ceq] = nonlinear(~, parts, ctx)
            x = parts.x;
            u = parts.u;
            p = ctx.p;
            mh = ctx.m;

            xhat = zeros(p + 1, ctx.nx);
            xhat(1, :) = ctx.x_init;
            uk = u(1, :);
            for i = 1 : p
                if i <= mh
                    uk = u(i, :);
                end
                xhat(i + 1, :) = ctx.model.step(x(i, :), uk);
            end

            sx = ctx.scaling.get('x', ctx.nx);
            ceq_blk = [(x(2:end, :) - xhat(2:end, :)) ./ sx;
                       (x(1, :)     - ctx.x_init)     ./ sx];
            ceq = reshape(ceq_blk, [], 1);
            c   = [];
        end
    end
end
