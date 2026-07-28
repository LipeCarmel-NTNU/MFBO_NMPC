function NMPC = build_nmpc(base, cfg, opts)
%BUILD_NMPC Construct and configure the NMPC controller for one theta.
%
%   NMPC = build_nmpc(base, cfg) builds an NMPC_terminal from the model and
%   dimensions in base and applies the horizons and weights in cfg (see
%   decode_theta). Setpoints and the terminal cost come from base unless
%   overridden.
%
%   Name-value options:
%     terminal_cost   'lqr'  terminal matrix from construct_P and base.LQR_data
%                     'zero' P = 0, the benchmark reference definition
%                     'none' leave P unset, for callers that assign it per
%                            setpoint segment
%     set_setpoint    copy base.xsp/base.usp into the controller (default true)
%     max_iter        fmincon MaxIterations (default base.optimizer_max_iter)
%     display         fmincon Display (default 'final-detailed')
%
%   FiniteDifferenceType is 'central' here as well as in the NMPC_terminal
%   class default. fmincon builds gradients by finite differences and forward
%   differences are too inaccurate in the flat regions of this problem, so the
%   setting is repeated at the point of use to keep it visible.

    arguments
        base struct
        cfg struct
        opts.terminal_cost (1,1) string {mustBeMember(opts.terminal_cost, ["lqr" "zero" "none"])} = "lqr"
        opts.set_setpoint (1,1) logical = true
        opts.max_iter double = []
        opts.display (1,1) string = "final-detailed"
    end

    if isempty(opts.max_iter)
        opts.max_iter = base.optimizer_max_iter;
    end

    NMPC = NMPC_terminal(base.model, base.nx, base.nu);

    %% Optimiser options
    % UseParallel only helps when a pool is already open; fmincon would
    % otherwise start one per solve.
    pool_empty = isempty(gcp('nocreate'));
    NMPC.optimizer_options.MaxIterations = opts.max_iter;
    NMPC.optimizer_options.UseParallel = ~pool_empty;
    NMPC.optimizer_options.Display = opts.display;
    NMPC.optimizer_options.FiniteDifferenceType = 'central';

    %% Horizons and weights
    NMPC.Ts = base.dt;
    NMPC.p = cfg.p;
    NMPC.m = cfg.m;
    NMPC.Q = cfg.Q;
    NMPC.Ru = cfg.Ru;
    NMPC.Rdu = cfg.Rdu;

    % Rebuild the decision-variable bounds for the new horizons.
    NMPC.constraints();

    %% Setpoints
    if opts.set_setpoint
        if isempty(base.xsp) || isempty(base.usp)
            error("base.xsp and base.usp are empty; build nmpc_base with set_setpoint true or pass set_setpoint false here.");
        end
        NMPC.xsp = base.xsp;
        NMPC.usp = base.usp;
    end

    %% Terminal cost
    switch opts.terminal_cost
        case "lqr"
            if isempty(base.LQR_data)
                error("base.LQR_data is empty; build nmpc_base with load_lqr true to use terminal_cost 'lqr'.");
            end
            NMPC.P = construct_P(base.LQR_data, NMPC.Q, NMPC.Ru, NMPC.Rdu);
        case "zero"
            NMPC.P = zeros(base.nx + base.nu);
        case "none"
            % The caller assigns P for each setpoint segment.
    end
end
