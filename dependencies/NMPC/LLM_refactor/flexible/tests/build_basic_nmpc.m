function nmpc = build_basic_nmpc(args)
    % build_basic_nmpc  Shared factory for the flexible-version tests.
    %
    %   Mirrors minimal/tests/build_basic_nmpc. Uses the same 2D decoupled
    %   stable LTI plant so test results can be cross-checked.

    arguments
        args.p         = 10
        args.m         = 5
        args.Ts        = 0.1
        args.y_sp      = [1 2]
        args.u_sp      = [0 0]
        args.Q         = diag([10 10])
        args.R         = diag([0.1 0.1])
        args.Ymin      = [-10 -10]
        args.Ymax      = [ 10  10]
        args.umin      = [-5 -5]
        args.umax      = [ 5  5]
        args.x_scale   = []
        args.u_scale   = []
        args.P         = []
        args.x_sp_term = []
        args.soft_mask = []
        args.rho_L1    = 0
        args.rho_L2    = 0
        args.S         = []
        args.dumax     = []
        args.log_enabled = false
    end

    nx = 2;
    nu = 2;
    f  = @(x, u) ((diag([-1 -2]) * x(:) + eye(2) * u(:)).');

    builder = NMPCBuilder() ...
        .withModel(ContinuousModel(f, nx, nu, RK4(args.Ts))) ...
        .withHorizon(args.p, args.m);

    if ~isempty(args.x_scale); builder.withScaling('x', args.x_scale); end
    if ~isempty(args.u_scale); builder.withScaling('u', args.u_scale); end
    if ~isempty(args.soft_mask)
        builder.withScaling('slack_state', ones(1, sum(logical(args.soft_mask))));
    end

    builder ...
        .addConstraint(StateBounds(args.Ymin, args.Ymax, soft=args.soft_mask)) ...
        .addConstraint(InputBounds(args.umin, args.umax));

    if ~isempty(args.dumax)
        builder.addConstraint(InputMovement(args.dumax));
    end

    builder.addCost(QuadraticTracking(args.Q, args.R, args.y_sp, args.u_sp));

    if ~isempty(args.P)
        xspT = args.x_sp_term;
        if isempty(xspT); xspT = args.y_sp; end
        builder.addCost(TerminalCost(args.P, xspT));
    end
    if ~isempty(args.S)
        builder.addCost(InputMovementCost(args.S));
    end
    if ~isempty(args.soft_mask) && (any(args.rho_L1 ~= 0) || any(args.rho_L2 ~= 0))
        builder.addCost(SlackPenalty(rho_L1=args.rho_L1, rho_L2=args.rho_L2));
    end

    if args.log_enabled
        builder.enableLog();
    end

    nmpc = builder.build();
end
