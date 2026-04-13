classdef NMPC < handle
    % NMPC  Composition-based nonlinear MPC orchestrator.
    %
    %   Built by NMPCBuilder. Wires a model, a frozen DecisionLayout, a
    %   Scaling object, an ordered list of constraints and costs, a solver
    %   backend, and (optionally) a HistoryLogger.
    %
    %   Contracts are duck-typed: constraints may implement any subset of
    %       register, linear_static, linear_per_call, nonlinear, bound_overrides
    %   Costs may implement
    %       register, evaluate
    %
    %   See LLM_refactor/flexible/STRUCTURE.md for the full design.

    properties (SetAccess = private)
        model
        p
        m
        nx
        nu
        ny

        layout          % DecisionLayout (frozen)
        scaling         % Scaling
        constraints     % cell array
        costs           % cell array
        backend         % solver backend

        log_enabled     % logical
        logger          % HistoryLogger (may be empty)

        % Cached pieces
        col_scale       % column scale factor: total×1 (w_phys = col_scale .* w_scaled)
        A_static        % static linear rows in scaled coordinates
        b_static
        wL_static       % scaled lower bound
        wU_static       % scaled upper bound

        % Runtime
        latest_wopt = []   % scaled
        latest_flag = NaN
    end

    methods
        function obj = NMPC(cfg)
            % Package `cfg` contains: model, p, m, ny, layout, scaling,
            % constraints (cell), costs (cell), backend, log_enabled, logger.
            obj.model       = cfg.model;
            obj.p           = cfg.p;
            obj.m           = cfg.m;
            obj.nx          = cfg.model.nx;
            obj.nu          = cfg.model.nu;
            obj.ny          = cfg.ny;
            obj.layout      = cfg.layout;
            obj.scaling     = cfg.scaling;
            obj.constraints = cfg.constraints;
            obj.costs       = cfg.costs;
            obj.backend     = cfg.backend;
            obj.log_enabled = cfg.log_enabled;
            obj.logger      = cfg.logger;

            obj.build_static_cache();
        end

        function log = get_log(obj)
            if isempty(obj.logger)
                log = struct([]);
            else
                log = obj.logger.records;
            end
        end

        function [uk, x_phys, u_phys, info] = solve(obj, x_init, u_prev)
            arguments
                obj
                x_init (1,:) double {mustBeFinite}
                u_prev (1,:) double {mustBeFinite}
            end
            assert(numel(x_init) == obj.nx, 'NMPC:x_init', ...
                'x_init must have length nx (%d).', obj.nx);
            assert(numel(u_prev) == obj.nu, 'NMPC:u_prev', ...
                'u_prev must have length nu (%d).', obj.nu);

            failed_before = (obj.latest_flag < 0);

            % 1. Cold / warm initial guess.
            if isempty(obj.latest_wopt) || failed_before
                w0 = obj.guess_from_initial(x_init, u_prev);
            else
                parts_prev = obj.layout.unpack(obj.latest_wopt, obj.scaling);
                u_prev_traj = parts_prev.u;
                u_shift = [u_prev_traj(2:end, :); u_prev_traj(end, :)];
                w0 = obj.guess_from_initial(x_init, u_shift);
            end

            [uk, x_phys, u_phys, fval] = obj.solve_once(w0, x_init, u_prev);

            % 2. Continuity-restoring retry on flag == -2.
            if obj.latest_flag == -2
                w0 = obj.guess_from_initial(x_init, u_phys);
                [uk, x_phys, u_phys, fval] = obj.solve_once(w0, x_init, u_prev);
            end

            % 3. Fallback.
            if obj.latest_flag < 0
                if failed_before
                    warning('NMPC:failed_again', 'NMPC failed twice in a row.');
                    uk = zeros(1, obj.nu);
                else
                    warning('NMPC:fallback', ...
                        'NMPC infeasible; falling back to previous trajectory.');
                    if isempty(obj.latest_wopt)
                        uk = zeros(1, obj.nu);
                    else
                        parts_prev = obj.layout.unpack(obj.latest_wopt, obj.scaling);
                        if obj.m > 1
                            uk = parts_prev.u(2, :);
                        else
                            uk = parts_prev.u(1, :);
                        end
                    end
                end
            end

            info = struct('flag', obj.latest_flag, 'fval', fval);
            if obj.log_enabled && ~isempty(obj.logger)
                rec = struct( ...
                    'x_init', x_init, 'u_prev', u_prev, 'uk', uk, ...
                    'w_opt',  obj.latest_wopt, 'flag', obj.latest_flag, ...
                    'fval',   fval, 't_wall', now);                       %#ok<TNOW1>
                obj.logger.append(rec);
            end
        end

        function x_next = step(obj, x, u)
            x_next = obj.model.step(x, u);
        end

        function w0 = guess_from_initial(obj, x_init, u_init)
            if size(u_init, 1) == 1
                u = [u_init; zeros(obj.m - 1, obj.nu)];
            else
                u = u_init;
            end
            u = max(u, 0);

            x = zeros(obj.p + 1, obj.nx);
            x(1, :) = x_init;
            uk = u(1, :);
            for i = 1 : obj.p
                if i <= obj.m
                    uk = u(i, :);
                end
                x(i + 1, :) = obj.model.step(x(i, :), uk);
            end

            parts = struct('x', x, 'u', u);
            % Zero-initialize any extra registered blocks (slacks, aux).
            names = obj.layout.list();
            for k = 1 : numel(names)
                nm = names{k};
                if ~isfield(parts, nm)
                    sz = obj.layout.size_of(nm);
                    parts.(nm) = zeros(sz);
                end
            end
            w0 = obj.layout.pack(parts, obj.scaling);
        end
    end

    methods (Access = private)
        function build_static_cache(obj)
            % Column scale factor: col_scale(j) such that w_phys_j = col_scale(j) * w_scaled_j.
            names = obj.layout.list();
            total = obj.layout.total_size();
            obj.col_scale = ones(total, 1);
            for k = 1 : numel(names)
                nm = names{k};
                sz = obj.layout.size_of(nm);     % [n_steps, n_cols]
                sc = obj.scaling.get(nm, sz(2)); % 1×n_cols
                blk = repmat(sc, sz(1), 1);      % n_steps×n_cols
                [lo, hi] = obj.layout.range_of(nm);
                obj.col_scale(lo : hi) = blk(:);
            end

            % Aggregate bound overrides (physical), then scale.
            lo_phys =  -inf(total, 1);
            hi_phys =   inf(total, 1);
            for k = 1 : numel(obj.constraints)
                c = obj.constraints{k};
                if ismethod(c, 'bound_overrides')
                    ov = c.bound_overrides(obj.layout);
                    fn = fieldnames(ov);
                    for j = 1 : numel(fn)
                        nm = fn{j};
                        if ~obj.layout.has(nm)
                            continue
                        end
                        [alo, ahi] = obj.layout.range_of(nm);
                        lo_phys(alo:ahi) = max(lo_phys(alo:ahi), ov.(nm).lo(:));
                        hi_phys(alo:ahi) = min(hi_phys(alo:ahi), ov.(nm).hi(:));
                    end
                end
            end
            obj.wL_static = lo_phys ./ obj.col_scale;
            obj.wU_static = hi_phys ./ obj.col_scale;

            % Aggregate static linear rows (in physical-w coords), then rescale columns.
            A_phys = zeros(0, total);
            b_phys = zeros(0, 1);
            for k = 1 : numel(obj.constraints)
                c = obj.constraints{k};
                if ismethod(c, 'linear_static')
                    [Ak, bk] = c.linear_static(obj.layout);
                    if ~isempty(Ak)
                        A_phys = [A_phys; Ak];       %#ok<AGROW>
                        b_phys = [b_phys; bk];       %#ok<AGROW>
                    end
                end
            end
            obj.A_static = A_phys .* obj.col_scale.';   % columns j scaled by col_scale(j)
            obj.b_static = b_phys;
        end

        function [uk, x_phys, u_phys, fval] = solve_once(obj, w0, x_init, u_prev)
            ctx = obj.make_ctx(x_init, u_prev);

            % Per-call linear rows.
            A = obj.A_static;
            b = obj.b_static;
            for k = 1 : numel(obj.constraints)
                c = obj.constraints{k};
                if ismethod(c, 'linear_per_call')
                    [Ak, bk] = c.linear_per_call(obj.layout, ctx);
                    if ~isempty(Ak)
                        Ak_scaled = Ak .* obj.col_scale.';
                        A = [A; Ak_scaled];     %#ok<AGROW>
                        b = [b; bk];            %#ok<AGROW>
                    end
                end
            end

            problem = struct( ...
                'objfun',  @(ws) obj.eval_cost(ws, ctx), ...
                'w0',      w0, ...
                'wL',      obj.wL_static, ...
                'wU',      obj.wU_static, ...
                'A',       A, 'b', b, ...
                'Aeq',     [], 'beq', [], ...
                'nonlcon', @(ws) obj.eval_nonlcon(ws, ctx), ...
                'options', []);

            [wopt_s, fval, flag] = obj.backend.solve(problem);

            if flag >= 0
                obj.latest_wopt = wopt_s;
            end
            obj.latest_flag = flag;

            parts  = obj.layout.unpack(wopt_s, obj.scaling);
            x_phys = parts.x;
            u_phys = parts.u;
            uk     = u_phys(1, :);
        end

        function J = eval_cost(obj, ws, ctx)
            parts = obj.layout.unpack(ws, obj.scaling);
            J = 0;
            for k = 1 : numel(obj.costs)
                c = obj.costs{k};
                if ismethod(c, 'evaluate')
                    J = J + c.evaluate(parts, ctx);
                end
            end
        end

        function [c_out, ceq_out] = eval_nonlcon(obj, ws, ctx)
            parts = obj.layout.unpack(ws, obj.scaling);
            c_out   = [];
            ceq_out = [];
            for k = 1 : numel(obj.constraints)
                c = obj.constraints{k};
                if ismethod(c, 'nonlinear')
                    [ck, ceqk] = c.nonlinear(parts, ctx);
                    c_out   = [c_out;   ck(:)];    %#ok<AGROW>
                    ceq_out = [ceq_out; ceqk(:)];  %#ok<AGROW>
                end
            end
        end

        function ctx = make_ctx(obj, x_init, u_prev)
            ctx = struct( ...
                'x_init',  x_init, ...
                'u_prev',  u_prev, ...
                'model',   obj.model, ...
                'scaling', obj.scaling, ...
                'layout',  obj.layout, ...
                'p',       obj.p, ...
                'm',       obj.m, ...
                'nx',      obj.nx, ...
                'nu',      obj.nu, ...
                'ny',      obj.ny);
        end
    end
end
