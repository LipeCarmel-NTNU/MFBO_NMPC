classdef NMPC < handle
    % NMPC  Single-class nonlinear MPC controller.
    %
    %   Fixed assumptions:
    %     - Continuous-time model xdot = f(x, u), single-shooting per Ts.
    %     - Multiple-shooting: each predicted state is a decision variable;
    %       continuity is enforced as nonlinear equality constraints.
    %     - Decision variables live in scaled coordinates (x_scale = 1 by
    %       default ⇒ identity, no behavioural difference).
    %     - fmincon as solver, RK4 as integrator.
    %
    %   Optional features (toggled by leaving the related property empty
    %   or zero ⇒ no contribution):
    %     - Terminal cost   .... P, x_sp_terminal
    %     - Soft state bnds .... soft_mask, rho_L1, rho_L2
    %     - Δu penalty ......... S
    %     - Δu hard bound ...... dumax
    %     - History log ........ log_enabled
    %
    %   Construction is by name-value pairs; see the constructor.
    %
    %   See LLM_refactor/minimal/STRUCTURE.md for the design rationale.

    properties
        %% Model & horizon (required)
        f                       % @(x,u) -> xdot   (1×nx out, x,u row vectors)
        nx                      % # states
        nu                      % # inputs
        ny                      % # outputs (defaults to nx)
        Ts                      % sampling time
        p                       % prediction horizon (steps)
        m                       % control horizon (steps), m <= p

        %% Tracking targets and weights (required)
        y_sp                    % 1×ny or p×ny setpoint
        u_sp                    % 1×nu or m×nu input reference
        Q                       % ny×ny output weight
        R                       % nu×nu input weight
        h_y = @(x) x            % output map; default y = x  (requires ny=nx)

        %% Bounds (required)
        Ymin                    % 1×nx
        Ymax                    % 1×nx
        umin                    % 1×nu
        umax                    % 1×nu

        %% Optional terminal cost
        P              = []     % nx×nx, empty ⇒ skipped
        x_sp_terminal  = []     % 1×nx; defaults to last row of y_sp if empty AND P is set

        %% Scaling (always present, default = ones)
        x_scale = []            % 1×nx, filled to ones in init
        u_scale = []            % 1×nu, filled to ones in init

        %% Optional soft state bounds
        soft_mask = []          % logical 1×nx, empty ⇒ no slacks
        rho_L1    = 0           % scalar or 1×n_soft
        rho_L2    = 0           % scalar or 1×n_soft

        %% Optional input-movement cost / hard bound
        S      = []             % nu×nu
        dumax  = []             % 1×nu (|Δu| ≤ dumax)

        %% Solver
        optimizer_options = optimoptions('fmincon', ...
            'Display','off','Algorithm','sqp', ...
            'MaxFunEvals',Inf,'MaxIterations',1000, ...
            'StepTolerance',1e-9,'OptimalityTolerance',1e-6, ...
            'ScaleProblem',true);

        %% Optional history log
        log_enabled = false
        log         = struct([])

        %% Runtime state
        latest_wopt = []        % stored in scaled coordinates
        latest_flag = NaN
    end

    properties (SetAccess = private)
        % Read-only layout/bookkeeping populated by init().
        n_soft   = 0            % number of soft states
        soft_idx = []           % indices into 1:nx of soft states
        len_x    = 0            % (p+1)*nx
        len_u    = 0            % m*nu
        len_s    = 0            % (p+1)*n_soft
    end

    properties (Access = private)
        % Cached linear-constraint pieces (in scaled-decision-vector space).
        A_soft = []             % soft state-bound rows (full, not sparse)
        b_soft = []
    end

    methods
        %% Construction
        function obj = NMPC(opts)
            % NMPC  Construct via name=value pairs.
            %
            %   nmpc = NMPC(f=@model, nx=4, nu=3, Ts=1/60, p=60, m=6, ...
            %               y_sp=ysp, u_sp=usp, Q=Q, R=R, ...
            %               Ymin=Ymin, Ymax=Ymax, umin=umin, umax=umax, ...);
            %
            %   Optional: ny, h_y, P, x_sp_terminal, x_scale, u_scale,
            %             soft_mask, rho_L1, rho_L2, S, dumax,
            %             optimizer_options, log_enabled.

            arguments
                % Required: model & horizon
                opts.f                  function_handle
                opts.nx           (1,1) double {mustBeInteger, mustBePositive}
                opts.nu           (1,1) double {mustBeInteger, mustBePositive}
                opts.Ts           (1,1) double {mustBePositive, mustBeFinite}
                opts.p            (1,1) double {mustBeInteger, mustBePositive}
                opts.m            (1,1) double {mustBeInteger, mustBePositive}

                % Required: tracking
                opts.y_sp         double {mustBeNonempty}
                opts.u_sp         double {mustBeNonempty}
                opts.Q            double {mustBeNonempty}
                opts.R            double {mustBeNonempty}

                % Required: bounds
                opts.Ymin   (1,:) double {mustBeNonempty}
                opts.Ymax   (1,:) double {mustBeNonempty}
                opts.umin   (1,:) double {mustBeNonempty}
                opts.umax   (1,:) double {mustBeNonempty}

                % Optional
                opts.ny                  double = []
                opts.h_y                 function_handle = @(x) x
                opts.P                   double = []
                opts.x_sp_terminal       double = []
                opts.x_scale       (1,:) double = []
                opts.u_scale       (1,:) double = []
                opts.soft_mask           = []           % logical or numeric mask
                opts.rho_L1              double = 0
                opts.rho_L2              double = 0
                opts.S                   double = []
                opts.dumax         (1,:) double = []
                opts.optimizer_options          = []
                opts.log_enabled   (1,1) logical = false
            end

            fn = fieldnames(opts);
            for k = 1 : numel(fn)
                v = opts.(fn{k});
                if strcmp(fn{k}, 'optimizer_options') && isempty(v)
                    continue        % keep class default
                end
                obj.(fn{k}) = v;
            end
            obj.init();
        end

        function init(obj)
            % Cross-field validation, default fill, cache static pieces.

            % Required fields. Name-value `arguments` entries are always
            % optional in MATLAB even without a default, so we re-check here.
            req = {'f','nx','nu','Ts','p','m','y_sp','u_sp','Q','R', ...
                   'Ymin','Ymax','umin','umax'};
            for k = 1 : numel(req)
                if isempty(obj.(req{k}))
                    error('NMPC:missing', ...
                        'Required property "%s" is empty.', req{k});
                end
            end

            if isempty(obj.ny);      obj.ny      = obj.nx;          end
            if isempty(obj.x_scale); obj.x_scale = ones(1, obj.nx); end
            if isempty(obj.u_scale); obj.u_scale = ones(1, obj.nu); end

            assert(numel(obj.x_scale) == obj.nx && all(obj.x_scale > 0), ...
                'NMPC:x_scale', 'x_scale must be 1×nx and strictly positive.');
            assert(numel(obj.u_scale) == obj.nu && all(obj.u_scale > 0), ...
                'NMPC:u_scale', 'u_scale must be 1×nu and strictly positive.');
            assert(obj.m <= obj.p, 'NMPC:horizon', 'Require m <= p.');
            assert(numel(obj.Ymin) == obj.nx && numel(obj.Ymax) == obj.nx, ...
                'NMPC:Ybounds', 'Ymin/Ymax must have length nx.');
            assert(numel(obj.umin) == obj.nu && numel(obj.umax) == obj.nu, ...
                'NMPC:ubounds', 'umin/umax must have length nu.');
            assert(isequal(size(obj.Q), [obj.ny obj.ny]), ...
                'NMPC:Qsize', 'Q must be ny×ny.');
            assert(isequal(size(obj.R), [obj.nu obj.nu]), ...
                'NMPC:Rsize', 'R must be nu×nu.');
            if ~isempty(obj.P)
                assert(isequal(size(obj.P), [obj.nx obj.nx]), ...
                    'NMPC:Psize', 'P must be nx×nx.');
            end
            if ~isempty(obj.S)
                assert(isequal(size(obj.S), [obj.nu obj.nu]), ...
                    'NMPC:Ssize', 'S must be nu×nu.');
            end
            if ~isempty(obj.dumax)
                assert(numel(obj.dumax) == obj.nu && all(obj.dumax > 0), ...
                    'NMPC:dumax', 'dumax must be 1×nu and strictly positive.');
            end

            % Slack bookkeeping
            if ~isempty(obj.soft_mask)
                obj.soft_mask = logical(obj.soft_mask(:).');
                assert(numel(obj.soft_mask) == obj.nx, ...
                    'soft_mask must have length nx.');
                obj.soft_idx = find(obj.soft_mask);
                obj.n_soft   = numel(obj.soft_idx);
            else
                obj.soft_idx = [];
                obj.n_soft   = 0;
            end

            obj.len_x = (obj.p + 1) * obj.nx;
            obj.len_u =  obj.m      * obj.nu;
            obj.len_s = (obj.p + 1) * obj.n_soft;

            % Build cached soft-bound linear rows (scaled space).
            obj.build_soft_rows();
        end

        %% Public solve  (mirrors NMPC_abstract.solve fallback flow)
        function [uk, x_phys, u_phys, info] = solve(obj, x_init, u_init)
            % SOLVE  Compute the next control move.
            %
            %   x_init : 1×nx current state (physical units)
            %   u_init : 1×nu previous applied input (physical units)

            arguments
                obj
                x_init (1,:) double {mustBeFinite}
                u_init (1,:) double {mustBeFinite}
            end
            assert(numel(x_init) == obj.nx, 'NMPC:x_init', ...
                'x_init must have length nx (%d).', obj.nx);
            assert(numel(u_init) == obj.nu, 'NMPC:u_init', ...
                'u_init must have length nu (%d).', obj.nu);

            failed_before = (obj.latest_flag < 0);  % NaN on startup -> false

            % 1. Cold solve or warm start
            if isempty(obj.latest_wopt) || failed_before
                w0 = obj.guess_from_initial(x_init, u_init);
                [uk, x_phys, u_phys, fval] = ...
                    obj.solve_optimization(w0, x_init, u_init);
            else
                % Shift previous u, keep last row, rebuild w0 by simulation.
                [~, u_prev] = obj.unpack_phys(obj.latest_wopt);
                u_shift = [u_prev(2:end, :); u_prev(end, :)];
                w0 = obj.guess_from_initial(x_init, u_shift);
                [uk, x_phys, u_phys, fval] = ...
                    obj.solve_optimization(w0, x_init, u_init);
            end

            % 2. Continuity-restoring retry on flag == -2
            if obj.latest_flag == -2
                w0 = obj.guess_from_initial(x_init, u_phys);
                [uk, x_phys, u_phys, fval] = ...
                    obj.solve_optimization(w0, x_init, u_init);
            end

            % 3. Fallback to u(k+1|k-1) or zeros
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
                        [~, u_prev] = obj.unpack_phys(obj.latest_wopt);
                        if obj.m > 1
                            uk = u_prev(2, :);
                        else
                            uk = u_prev(1, :);
                        end
                    end
                end
            end

            % 4. Log
            info = struct('flag', obj.latest_flag, 'fval', fval);
            if obj.log_enabled
                obj.append_log(x_init, u_init, uk, fval);
            end
        end

        %% Internals
        function [uk, x_phys, u_phys, fval] = solve_optimization(obj, w0, x_init, u_init)
            % w0 is already in scaled-decision-vector form (pack_phys
            % converts physical inputs to scaled w internally).

            % Bounds (scaled).
            [wL, wU] = obj.bounds_scaled();

            % Linear inequality rows (scaled): static soft + per-call Δu.
            [A, b] = obj.linear_ineq(u_init);

            [wopt_s, fval, exitflag] = fmincon( ...
                @(ws) obj.objfun(ws, u_init), w0, ...
                A, b, [], [], wL, wU, ...
                @(ws) obj.confun(ws, x_init), ...
                obj.optimizer_options);

            if exitflag >= 0
                obj.latest_wopt = wopt_s;        % store scaled
            end
            obj.latest_flag = exitflag;

            [x_phys, u_phys] = obj.unpack_phys(wopt_s);
            uk = u_phys(1, :);
        end

        %% Objective (scaled inputs)
        function J = objfun(obj, ws, u_prev)
            [x, u, s] = obj.unpack_phys(ws);

            J = 0;

            % Tracking on outputs over the prediction horizon
            ysp = obj.expand_setpoint(obj.y_sp, obj.p, obj.ny);
            for k = 1 : obj.p
                yk = obj.h_y(x(k, :));
                e  = (yk - ysp(k, :)).';
                J  = J + e.' * obj.Q * e;
            end

            % Input reference penalty over the control horizon
            usp = obj.expand_setpoint(obj.u_sp, obj.m, obj.nu);
            for k = 1 : obj.m
                ek = (u(k, :) - usp(k, :)).';
                J  = J + ek.' * obj.R * ek;
            end

            % Optional Δu cost
            if ~isempty(obj.S)
                du = diff([u_prev; u], [], 1);
                for k = 1 : obj.m
                    duk = du(k, :).';
                    J = J + duk.' * obj.S * duk;
                end
            end

            % Optional terminal cost
            if ~isempty(obj.P)
                if isempty(obj.x_sp_terminal)
                    xsp_term = ysp(end, :);          % assumes ny == nx
                else
                    xsp_term = obj.x_sp_terminal;
                end
                e = (x(end, :) - xsp_term).';
                J = J + e.' * obj.P * e;
            end

            % Optional slack penalty (L1 + L2)
            if ~isempty(s)
                [N_, ~] = size(s);
                sv = s(:);
                wL1 = NMPC.tile_weights(obj.rho_L1, N_, obj.n_soft);
                J = J + sum(wL1 .* sv);
                if any(obj.rho_L2(:) ~= 0)
                    wL2 = NMPC.tile_weights(obj.rho_L2, N_, obj.n_soft);
                    J = J + sv.' * (wL2 .* sv);
                end
            end
        end

        %% Continuity equality constraints
        function [c, ceq] = confun(obj, ws, x_init)
            [x, u, ~] = obj.unpack_phys(ws);

            xhat = zeros(obj.p + 1, obj.nx);
            xhat(1, :) = x_init;
            uk = u(1, :);
            for i = 1 : obj.p
                if i <= obj.m
                    uk = u(i, :);
                end
                xhat(i + 1, :) = obj.step(x(i, :), uk);
            end

            % Initial-state and continuity equalities, scaled by x_scale
            % so residuals are dimensionless (matches NMPC_abstract_scaled).
            sx = obj.x_scale;
            ceq_blk = [(x(2:end, :) - xhat(2:end, :)) ./ sx;
                       (x(1, :)     - x_init)         ./ sx];
            ceq = reshape(ceq_blk, [], 1);
            c   = [];
        end

        %% One sampling step (RK4)
        function x_next = step(obj, x, u)
            h = obj.Ts;
            k1 = obj.fr(x,            u);
            k2 = obj.fr(x + 0.5*h*k1, u);
            k3 = obj.fr(x + 0.5*h*k2, u);
            k4 = obj.fr(x +     h*k3, u);
            x_next = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        end

        function v = fr(obj, x, u)
            % Force the model output to a row vector so RK4 broadcasts cleanly
            % regardless of the user's column/row convention.
            v = reshape(obj.f(x, u), 1, []);
        end

        %% Initial guess (physical)
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
                x(i + 1, :) = obj.step(x(i, :), uk);
            end

            if obj.n_soft > 0
                s = zeros(obj.p + 1, obj.n_soft);
            else
                s = [];
            end
            w0 = obj.pack_phys(x, u, s);
        end

        %% Pack / unpack with scaling
        function w = pack_phys(obj, x, u, s)
            % Pack PHYSICAL (x, u, s) arrays into the SCALED decision vector
            % w used by fmincon and the internal cost/constraint functions.
            xs = x ./ obj.x_scale;
            us = u ./ obj.u_scale;
            if isempty(s)
                w = [reshape(xs, [], 1); reshape(us, [], 1)];
            else
                ss = s ./ obj.x_scale(obj.soft_idx);
                w  = [reshape(xs, [], 1); reshape(us, [], 1); reshape(ss, [], 1)];
            end
        end

        function [x, u, s] = unpack_phys(obj, ws)
            % Unpack the SCALED decision vector ws back into PHYSICAL arrays.
            xs = reshape(ws(1 : obj.len_x),                       [], obj.nx);
            us = reshape(ws(obj.len_x + (1 : obj.len_u)),         [], obj.nu);
            x  = xs .* obj.x_scale;
            u  = us .* obj.u_scale;
            if obj.len_s > 0
                ss = reshape(ws(obj.len_x + obj.len_u + (1 : obj.len_s)), ...
                             [], obj.n_soft);
                s  = ss .* obj.x_scale(obj.soft_idx);
            else
                s = [];
            end
        end

        %% Bounds (scaled)
        function [wL, wU] = bounds_scaled(obj)
            % State bounds: relax soft entries to ±Inf (linear rows handle them).
            xL = repmat(obj.Ymin, obj.p + 1, 1);
            xU = repmat(obj.Ymax, obj.p + 1, 1);
            if obj.n_soft > 0
                xL(:, obj.soft_idx) = -Inf;
                xU(:, obj.soft_idx) =  Inf;
            end
            xL = xL ./ obj.x_scale;
            xU = xU ./ obj.x_scale;

            uL = repmat(obj.umin, obj.m, 1) ./ obj.u_scale;
            uU = repmat(obj.umax, obj.m, 1) ./ obj.u_scale;

            wL = [reshape(xL, [], 1); reshape(uL, [], 1)];
            wU = [reshape(xU, [], 1); reshape(uU, [], 1)];

            if obj.len_s > 0
                wL = [wL; zeros(obj.len_s, 1)];     % s >= 0
                wU = [wU;  inf(obj.len_s, 1)];
            end
        end

        %% Linear inequalities (scaled)
        function [A, b] = linear_ineq(obj, u_prev)
            A = obj.A_soft;
            b = obj.b_soft;

            if ~isempty(obj.dumax)
                [Adu, bdu] = obj.build_du_rows(u_prev);
                A = [A; Adu];
                b = [b; bdu];
            end
        end

        %% Build cached soft-bound rows (scaled, sparse)
        function build_soft_rows(obj)
            if obj.n_soft == 0
                obj.A_soft = zeros(0, obj.len_x + obj.len_u + obj.len_s);
                obj.b_soft = zeros(0, 1);
                return
            end
            N    = obj.p + 1;
            ns   = obj.n_soft;
            cols = obj.len_x + obj.len_u + obj.len_s;

            % Indices into ws of x(:, soft_idx) and s(:, :).
            xs_lin = zeros(N, ns);
            for j = 1 : ns
                col = obj.soft_idx(j);
                xs_lin(:, j) = (col - 1) * N + (1 : N).';
            end
            ss_lin = obj.len_x + obj.len_u + reshape(1 : obj.len_s, N, ns);

            % Soft-state scale factors (one per soft state, broadcast over N).
            sxs = obj.x_scale(obj.soft_idx);

            rows_per_block = N * ns;

            % Upper:  x_phys - s_phys ≤ Ymax
            %  ⇒    sx*xs - sx*ss   ≤ Ymax
            r_idx  = (1 : rows_per_block).';
            sx_blk = repmat(sxs, N, 1);                  % N×ns column-major
            sx_blk = sx_blk(:);

            A_up = sparse([r_idx; r_idx], ...
                          [xs_lin(:); ss_lin(:)], ...
                          [sx_blk;   -sx_blk], ...
                          rows_per_block, cols);
            b_up = repmat(obj.Ymax(obj.soft_idx).', N, 1);

            % Lower: -x_phys - s_phys ≤ -Ymin
            A_lo = sparse([r_idx; r_idx], ...
                          [xs_lin(:); ss_lin(:)], ...
                          [-sx_blk;  -sx_blk], ...
                          rows_per_block, cols);
            b_lo = -repmat(obj.Ymin(obj.soft_idx).', N, 1);

            % fmincon's sqp algorithm doesn't support sparse A; convert.
            obj.A_soft = full([A_up; A_lo]);
            obj.b_soft = [b_up; b_lo];
        end

        %% Build per-call Δu rows (scaled)
        function [A, b] = build_du_rows(obj, u_prev)
            % |u(k) - u(k-1)| ≤ dumax, with u(0) = u_prev (physical).
            cols = obj.len_x + obj.len_u + obj.len_s;
            nu_  = obj.nu;
            mh   = obj.m;

            % Dense; matrix is small (2*m*nu × cols) and SQP requires it.
            A = zeros(2 * mh * nu_, cols);
            b = zeros(2 * mh * nu_, 1);

            row = 0;
            for k = 1 : mh
                for j = 1 : nu_
                    % u is stored column-major as m × nu, so the linear
                    % index of u(k,j) inside the u block is (j-1)*m + k.
                    col_k = obj.len_x + (j - 1) * mh + k;
                    sj    = obj.u_scale(j);

                    if k == 1
                        % u(1) - u_prev:  sj * us(1,j) - u_prev(j) ≤ dumax(j)
                        row = row + 1;
                        A(row, col_k) =  sj;
                        b(row)        =  obj.dumax(j) + u_prev(j);

                        row = row + 1;
                        A(row, col_k) = -sj;
                        b(row)        =  obj.dumax(j) - u_prev(j);
                    else
                        col_km1 = obj.len_x + (j - 1) * mh + (k - 1);
                        row = row + 1;
                        A(row, col_k)   =  sj;
                        A(row, col_km1) = -sj;
                        b(row)          =  obj.dumax(j);

                        row = row + 1;
                        A(row, col_k)   = -sj;
                        A(row, col_km1) =  sj;
                        b(row)          =  obj.dumax(j);
                    end
                end
            end
        end

        %% History log
        function append_log(obj, x_init, u_init, uk, fval)
            rec = struct( ...
                'x_init',     x_init, ...
                'u_init',     u_init, ...
                'uk',         uk, ...
                'y_sp',       obj.y_sp, ...
                'u_sp',       obj.u_sp, ...
                'x_sp_term',  obj.x_sp_terminal, ...
                'Q',          obj.Q, ...
                'R',          obj.R, ...
                'P',          obj.P, ...
                'S',          obj.S, ...
                'rho_L1',     obj.rho_L1, ...
                'rho_L2',     obj.rho_L2, ...
                'w_opt',      obj.latest_wopt, ...
                'flag',       obj.latest_flag, ...
                'fval',       fval, ...
                't_wall',     now);                         %#ok<TNOW1>
            if isempty(obj.log)
                obj.log = rec;
            else
                obj.log(end + 1) = rec;
            end
        end
    end

    methods (Static)
        function M = expand_setpoint(sp, n_rows, n_cols)
            % Accept 1×ncols (broadcast) or n_rows×ncols.
            if size(sp, 1) == 1
                M = repmat(sp, n_rows, 1);
            else
                assert(size(sp, 1) == n_rows && size(sp, 2) == n_cols, ...
                    'Setpoint has incompatible size.');
                M = sp;
            end
        end

        function w = tile_weights(rho, N, ns)
            % Expand a slack-penalty weight to a length-(N*ns) column vector.
            % Accepts: scalar, 1×ns (per soft state), or already (N*ns)×1.
            if isscalar(rho)
                w = rho * ones(N * ns, 1);
            elseif numel(rho) == ns
                w = reshape(repmat(rho(:).', N, 1), [], 1);
            elseif numel(rho) == N * ns
                w = rho(:);
            else
                error('NMPC:rho_size', ...
                      'Slack weight must be scalar, 1×n_soft, or (p+1)*n_soft.');
            end
        end
    end
end
