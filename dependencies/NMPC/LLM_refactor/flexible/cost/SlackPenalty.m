classdef SlackPenalty < handle
    % SlackPenalty  L1 (exact) plus optional L2 on the `slack_state` block.
    %
    %   c = SlackPenalty(rho_L1=ρ)
    %   c = SlackPenalty(rho_L1=ρ1, rho_L2=ρ2)
    %
    %   No-op if no slack block was registered.

    properties (SetAccess = private)
        rho_L1
        rho_L2
        block
    end

    methods
        function obj = SlackPenalty(args)
            arguments
                args.rho_L1 double = 0
                args.rho_L2 double = 0
                args.block         = 'slack_state'
            end
            obj.rho_L1 = args.rho_L1;
            obj.rho_L2 = args.rho_L2;
            obj.block  = char(args.block);
        end

        function J = evaluate(obj, parts, ~)
            if ~isfield(parts, obj.block)
                J = 0;
                return
            end
            s  = parts.(obj.block);
            sv = s(:);
            J  = 0;
            if any(obj.rho_L1(:) ~= 0)
                w1 = SlackPenalty.tile(obj.rho_L1, size(s));
                J  = J + sum(w1 .* sv);
            end
            if any(obj.rho_L2(:) ~= 0)
                w2 = SlackPenalty.tile(obj.rho_L2, size(s));
                J  = J + sv.' * (w2 .* sv);
            end
        end
    end

    methods (Static, Access = private)
        function w = tile(rho, sz)
            N  = sz(1);
            ns = sz(2);
            if isscalar(rho)
                w = rho * ones(N * ns, 1);
            elseif numel(rho) == ns
                w = reshape(repmat(rho(:).', N, 1), [], 1);
            elseif numel(rho) == N * ns
                w = rho(:);
            else
                error('SlackPenalty:rho_size', ...
                    'Weight must be scalar, 1×n_soft, or (p+1)·n_soft.');
            end
        end
    end
end
