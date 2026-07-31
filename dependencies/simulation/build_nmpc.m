function nmpc = build_nmpc(base, cfg, opts)
%BUILD_NMPC Construct and configure the NMPC controller for one theta.
%
%   nmpc = build_nmpc(base, cfg) constructs an NMPC from the model, dimensions,
%   bounds and scaling in base and the horizons and weights in cfg (see
%   decode_theta). Setpoints and the terminal cost come from base unless
%   overridden. This is the only place in the pipeline that calls the NMPC
%   constructor.
%
%   The controller takes every field at construction, so nothing here needs a
%   later rebuild of the decision-variable bounds. The state bounds are relaxed
%   by an L1 slack wherever base.soft_mask is true; the slack carries the linear
%   penalty rho_L1 and no quadratic term.
%
%   Name-value options:
%     terminal_cost   'lqr'  terminal matrix from construct_P and base.LQR_data
%                     'zero' P = 0, the benchmark reference definition
%                     'none' leave P unset, for callers that assign it per
%                            setpoint segment together with x_term and u_term
%     set_setpoint    copy base.xsp/base.usp into the controller (default true)
%     max_iter        fmincon MaxIterations (default base.optimizer_max_iter)
%     display         fmincon Display (default 'final-detailed')
%     rho_L1          L1 slack penalty (default base.rho_L1)
%
%   FiniteDifferenceType is 'central' because fmincon builds gradients by
%   finite differences and forward differences are too inaccurate in the flat
%   regions of this problem. StepTolerance stays at the class default of 1e-7.

    arguments
        base struct
        cfg struct
        opts.terminal_cost (1,1) string {mustBeMember(opts.terminal_cost, ["lqr" "zero" "none"])} = "lqr"
        opts.set_setpoint (1,1) logical = true
        opts.max_iter double = []
        opts.display (1,1) string = "final"
        opts.rho_L1 double = []
    end

    if isempty(opts.max_iter)
        opts.max_iter = base.optimizer_max_iter;
    end
    if isempty(opts.rho_L1)
        opts.rho_L1 = base.rho_L1;
    end

    %% Setpoints
    % The controller requires a setpoint at construction. A caller that passes
    % set_setpoint false recomputes it for every segment, so the placeholder
    % below is overwritten by the setpoint hook before the first solve. Pair
    % set_setpoint false with the setpoint_fn option of nmpc_run_case.
    if opts.set_setpoint
        if isempty(base.xsp) || isempty(base.usp)
            error("base.xsp and base.usp are empty; build nmpc_base with set_setpoint true or pass set_setpoint false here.");
        end
        x_sp = base.xsp;
        u_sp = base.usp;
    else
        x_sp = zeros(1, base.nx);
        u_sp = zeros(1, base.nu);
    end

    %% Model
    % The controller integrates xdot(x, u) with a row state and expects a
    % 1-by-nx return, while base.model follows the orientation of its input.
    xdot = @(x, u) reshape(base.model(x(:), u), 1, []);

    %% Construction
    args = { ...
        'xdot', xdot, 'nx', base.nx, 'nu', base.nu, 'Ts', base.dt, ...
        'p', cfg.p, 'm', cfg.m, ...
        'x_sp', x_sp, 'u_sp', u_sp, ...
        'Q', cfg.Q, 'R_u', cfg.Ru, 'R_du', cfg.Rdu, ...
        'Xmin', base.Xmin, 'Xmax', base.Xmax, ...
        'umin', base.umin, 'umax', base.umax, ...
        'x_scale', base.x_scale, 'u_scale', base.u_scale, ...
        'soft_mask', base.soft_mask, 'rho_L1', opts.rho_L1};

    %% Terminal cost
    % P is augmented, (nx+nu)-by-(nx+nu), on z = [x_end - x_term; u_last -
    % u_term]. The terminal reference is the tracking setpoint, which is the
    % point construct_P linearises about.
    switch opts.terminal_cost
        case "lqr"
            if isempty(base.LQR_data)
                error("base.LQR_data is empty; build nmpc_base with load_lqr true to use terminal_cost 'lqr'.");
            end
            P = construct_P(base.LQR_data, cfg.Q, cfg.Ru, cfg.Rdu);
        case "zero"
            P = zeros(base.nx + base.nu);
        case "none"
            % The caller assigns P, x_term and u_term for each segment.
            P = [];
    end
    if ~isempty(P)
        args = [args, {'P', P, 'x_term', x_sp, 'u_term', u_sp}];
    end

    nmpc = NMPC(args{:});

    %% Optimiser options
    % UseParallel only helps when a pool is already open; fmincon would
    % otherwise start one per solve.
    pool_empty = isempty(gcp('nocreate'));
    nmpc.optimizer_options.MaxIterations = opts.max_iter;
    nmpc.optimizer_options.UseParallel = ~pool_empty;
    nmpc.optimizer_options.Display = opts.display;
    nmpc.optimizer_options.FiniteDifferenceType = 'central';
end
