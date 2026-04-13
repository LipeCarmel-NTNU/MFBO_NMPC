classdef StateBounds < handle
    % StateBounds  Hard or per-state-soft bounds on the `x` block.
    %
    %   c = StateBounds(Ymin, Ymax)                    % all hard
    %   c = StateBounds(Ymin, Ymax, soft=mask)         % soft mask (1×nx logical)
    %
    %   When any state is soft, a slack block `slack_state` of shape
    %   (p+1)×n_soft is registered. Bound inequalities for soft states are
    %   implemented as linear rows (in scaled w):
    %       x_i - s_i ≤ Ymax_i ,  -x_i - s_i ≤ -Ymin_i ,  s_i ≥ 0.
    %
    %   Hard bounds become `bound_overrides` applied to wL/wU.

    properties (SetAccess = private)
        Ymin
        Ymax
        soft         % logical 1×nx (all false ⇒ hard everywhere)
        soft_idx
        n_soft
    end

    methods
        function obj = StateBounds(Ymin, Ymax, args)
            arguments
                Ymin (1,:) double
                Ymax (1,:) double
                args.soft = []
            end
            assert(numel(Ymin) == numel(Ymax), ...
                'StateBounds:size', 'Ymin/Ymax size mismatch.');
            obj.Ymin = Ymin;
            obj.Ymax = Ymax;
            if isempty(args.soft)
                obj.soft = false(1, numel(Ymin));
            else
                obj.soft = logical(args.soft(:).');
                assert(numel(obj.soft) == numel(Ymin), ...
                    'StateBounds:soft', 'soft mask size mismatch with Ymin.');
            end
            obj.soft_idx = find(obj.soft);
            obj.n_soft   = numel(obj.soft_idx);
        end

        function register(obj, layout)
            % Register the slack block if any states are soft. Needs `x`
            % block to already be registered so we can read p+1.
            if obj.n_soft == 0
                return
            end
            assert(layout.has('x'), 'StateBounds:x_missing', ...
                'Register the x block before StateBounds.');
            sz = layout.size_of('x');    % [p+1, nx]
            layout.add('slack_state', sz(1), obj.n_soft);
        end

        function ov = bound_overrides(obj, layout)
            % Hard states -> tight wL/wU. Soft states -> ±Inf (rows handle).
            % Slacks -> [0, Inf).
            sz   = layout.size_of('x');  % [p+1, nx]
            nx   = sz(2);
            N    = sz(1);

            lo = repmat(obj.Ymin, N, 1);
            hi = repmat(obj.Ymax, N, 1);
            if obj.n_soft > 0
                lo(:, obj.soft_idx) = -Inf;
                hi(:, obj.soft_idx) =  Inf;
            end

            ov.x = struct('lo', lo, 'hi', hi);
            if obj.n_soft > 0
                ov.slack_state = struct( ...
                    'lo', zeros(N, obj.n_soft), ...
                    'hi',   inf(N, obj.n_soft));
            end
        end

        function [A, b] = linear_static(obj, layout)
            % Scaled rows: sx*xs - sx_s*ss ≤ Ymax   (and lower)
            if obj.n_soft == 0
                A = zeros(0, layout.total_size());
                b = zeros(0, 1);
                return
            end
            cols = layout.total_size();

            x_lin = layout.linear_indices('x');           % (p+1)×nx
            s_lin = layout.linear_indices('slack_state'); % (p+1)×n_soft
            N     = size(x_lin, 1);

            % Scale for x and slack columns comes via ctx-less build; use
            % ctx in linear_per_call if scaling needs to be applied. For
            % static rows we assume the layout is in physical pack units —
            % scaling is injected through ctx in the orchestrator by
            % substituting per-column scales. To keep this self-contained,
            % we emit rows in physical coordinates and let the orchestrator
            % rescale columns with scaling.get().
            %
            % But here we need the rows in scaled-w coordinates at the end.
            % Delegate the column-scaling multiplication to the orchestrator
            % via a convention: linear_static returns rows assuming identity
            % scaling; the orchestrator rescales each column later.

            rows_per = N * obj.n_soft;
            r = (1 : rows_per).';

            xs_col = x_lin(:, obj.soft_idx);    % N×n_soft
            ss_col = s_lin;                      % N×n_soft

            % Upper:  xs - ss ≤ Ymax(soft)
            A_up = sparse([r; r], ...
                          [xs_col(:); ss_col(:)], ...
                          [ ones(rows_per,1); -ones(rows_per,1)], ...
                          rows_per, cols);
            b_up = repmat(obj.Ymax(obj.soft_idx).', N, 1);

            % Lower: -xs - ss ≤ -Ymin(soft)
            A_lo = sparse([r; r], ...
                          [xs_col(:); ss_col(:)], ...
                          [-ones(rows_per,1); -ones(rows_per,1)], ...
                          rows_per, cols);
            b_lo = -repmat(obj.Ymin(obj.soft_idx).', N, 1);

            A = full([A_up; A_lo]);
            b = [b_up; b_lo];
        end
    end
end
