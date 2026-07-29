function base = nmpc_base(opts)
%NMPC_BASE One-time initialisation shared across all theta evaluations.
%
%   base = nmpc_base() builds the struct that every simulation entry point
%   reuses: model and plant handles, dimensions, the sampling grid, a fixed
%   measurement-noise realisation, and, when requested, the steady-state
%   setpoints and the LQR data needed to construct the terminal cost.
%
%   Name-value options:
%     sigma_y           measurement noise std dev per state (default
%                       [0.001 0.1 0.1]; pass [0 0 0] for a noise-free run)
%     Ts                sampling time in hours (default 1/60)
%     tf                horizon in hours (default 10)
%     Vsp, Xsp          setpoints passed to find_ss (default 1 and 20)
%     set_setpoint      store xsp/usp in base (default true). The
%                       setpoint-schedule run sets this false because it
%                       recomputes the setpoint for every segment.
%     load_lqr          load LQR_data.mat for construct_P (default true)
%     phi_coeffs_path   when non-empty, load the fitted fidelity surrogate
%                       coefficients from this file into base.phi
%     optimizer_max_iter  fmincon MaxIterations (default 100)
%
%   The noise realisation depends on the random state at call time, so seed
%   with rng before calling to keep runs comparable.

    arguments
        opts.sigma_y (1,:) double = [0.001 0.1 0.1]
        opts.Ts (1,1) double {mustBePositive} = 1/60
        opts.tf (1,1) double {mustBePositive} = 10
        opts.Vsp (1,1) double = 1
        opts.Xsp (1,1) double = 20
        opts.set_setpoint (1,1) logical = true
        opts.load_lqr (1,1) logical = true
        opts.phi_coeffs_path string = ""
        opts.optimizer_max_iter (1,1) double = 100
    end

    %% Project layout
    % This file lives in <root>/dependencies/simulation, so the project root
    % is two levels up. Resolving it here keeps the initialiser independent of
    % the current working directory.
    this_dir = fileparts(mfilename('fullpath'));
    base.project_root = fileparts(fileparts(this_dir));

    %% Model and plant
    % Non-negativity is enforced for states 2 and 3 during ODE integration.
    base.ode_opt = odeset('NonNegative', [2 3]);

    % get_par is a script that defines par in this workspace.
    get_par
    base.par = par;

    base.model = @(x, u) dilution_reduced(0, x, u(:)', base.par);
    base.plant = base.model;

    base.nx = 3;
    base.nu = 3;

    %% Sampling grid
    base.dt = opts.Ts;
    base.tf = opts.tf;
    base.N = ceil(base.tf / base.dt) + 1;
    base.T = (0:base.N-1).' * base.dt;
    base.tspan = [0 base.dt];               % one-step integration window

    %% Measurement noise
    % One realisation is drawn here and reused by every theta evaluation, so
    % controllers are compared against the same disturbance sequence.
    base.sigma_y = opts.sigma_y;
    base.noise = randn(base.N, base.nx) .* base.sigma_y;

    %% Setpoints
    base.V_sp = opts.Vsp;
    base.X_sp = opts.Xsp;
    if opts.set_setpoint
        [xss, uss] = find_ss(base.V_sp, base.X_sp, base.par, base.model, base.ode_opt);
        base.xsp = xss;
        base.usp = uss;
    else
        base.xsp = [];
        base.usp = [];
    end

    %% Terminal cost data
    if opts.load_lqr
        lqr_path = fullfile(base.project_root, 'LQR_data.mat');
        if ~isfile(lqr_path)
            error("LQR_data.mat not found at %s. Run TerminalCost/main_TerminalP.m first.", lqr_path);
        end
        S = load(lqr_path, 'LQR_data');
        base.LQR_data = S.LQR_data;
    else
        base.LQR_data = [];
    end

    %% Optimiser budget
    base.optimizer_max_iter = opts.optimizer_max_iter;

    %% Fidelity surrogate
    % Shape parameters of the cost fraction phi(z) = I_z(a, b), one pair for
    % each objective. The fit uses the runs collected so far. This field stays
    % empty during initialization, where the driver reports the cost measured at
    % the simulated fidelity and applies no scaling.
    %
    % The BO loop refits phi during the run. main_BO.m reloads this field before
    % every evaluation, so base.phi.vintage names the fit that scaled a row.
    if strlength(opts.phi_coeffs_path) > 0
        base.phi = load_phi_coeffs(opts.phi_coeffs_path);
    else
        base.phi = [];
    end
end
