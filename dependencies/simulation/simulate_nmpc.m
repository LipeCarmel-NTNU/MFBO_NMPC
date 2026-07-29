function out = simulate_nmpc(base, theta, opts)
%SIMULATE_NMPC Evaluate one theta over a set of initial-condition cases.
%
%   out = simulate_nmpc(base, theta) decodes theta, builds the controller,
%   runs every case, and aggregates the costs and runtimes.
%
%   Horizon
%     opts.horizon = "fidelity" runs tf = base.tf * f, where f is theta(1).
%       This is the fidelity used by the optimisation loop: a low f buys a
%       cheap, truncated evaluation.
%     opts.horizon = "full" ignores f and runs the whole base.tf, which is
%       what the benchmark and the setpoint-schedule comparison need.
%
%   Scaling by phi
%     opts.extrapolate divides the measured partial costs by the fitted cost
%     fraction phi(z) = I_z(a, b). This estimates the full-horizon totals. It
%     needs base.phi. It stays off during initialization, where out.SSE and
%     out.SSdU are the costs measured at the simulated fidelity. The fit of phi
%     uses those measured costs.
%
%     The measured costs survive the scaling as out.SSE_measured and
%     out.SSdU_measured. The function also stores the two divisors and the
%     surrogate vintage. You can therefore recompute the full-horizon estimate
%     under a later fit without a new simulation.
%
%   Timing
%     out.wall_s holds the wall time of each stage of this function.
%
%   Name-value options:
%     horizon        "fidelity" (default) or "full"
%     extrapolate    estimate full-horizon costs from base.phi (default false)
%     terminal_cost  "lqr" (default), "zero" or "none"; see build_nmpc
%     set_setpoint   copy base.xsp/base.usp into the controller (default true)
%     x0             one row per case (default [1.0 10 0; 1.1 25 5])
%     setpoint_fn    per-step setpoint hook; see nmpc_run_case
%     verbosity      "full" (default), "control" or "quiet"
%     run_id         identifier for the log and for checkpoint validation
%     log_path       SIMULATIONS_LOG.txt path; empty disables flag logging
%     checkpoint_path  partial .mat path; empty disables checkpoint and resume
%     max_iter       fmincon MaxIterations (default base.optimizer_max_iter)

    arguments
        base struct
        theta (1,:) double
        opts.horizon (1,1) string {mustBeMember(opts.horizon, ["fidelity" "full"])} = "fidelity"
        opts.extrapolate (1,1) logical = false
        opts.terminal_cost (1,1) string {mustBeMember(opts.terminal_cost, ["lqr" "zero" "none"])} = "lqr"
        opts.set_setpoint (1,1) logical = true
        opts.x0 (:,:) double = [1.0, 10, 0; 1.1, 25, 5]
        opts.setpoint_fn = []
        opts.verbosity (1,1) string {mustBeMember(opts.verbosity, ["full" "control" "quiet"])} = "full"
        opts.run_id (1,1) string = ""
        opts.log_path (1,1) string = ""
        opts.checkpoint_path (1,1) string = ""
        opts.max_iter double = []
    end

    t_call = tic;
    cfg = decode_theta(theta, base.nx, base.nu);
    if opts.verbosity ~= "quiet"
        disp('Run cfg:')
        disp(cfg)
    end

    %% Horizon
    if opts.horizon == "fidelity"
        tf = base.tf * cfg.f;
    else
        tf = base.tf;
    end
    N = ceil(tf / base.dt) + 1;
    N = min(N, base.N);

    %% Fidelity surrogate
    % phi(z) is the share of the full-horizon cost that the run accumulates by
    % fidelity z. I_z(a, b) already lies in [0, 1] and equals 1 at z = 1, so the
    % floor is the only guard that the caller needs. The floor caps the scaling
    % at 100x, so a very small phi cannot turn a small partial cost into a huge
    % estimate. Every evaluation records whether the floor acted.
    t_phi = tic;
    phi_floor = 0.01;
    phi_vintage = NaN;
    phi_floored = false;

    if opts.extrapolate
        if isempty(base.phi)
            error(['opts.extrapolate is true but base.phi is empty. Build ' ...
                   'nmpc_base with phi_coeffs_path set, or reload the field ' ...
                   'with load_phi_coeffs.']);
        end
        raw_SSE = phi_eval(cfg.f, base.phi.SSE.a, base.phi.SSE.b);
        raw_SSdU = phi_eval(cfg.f, base.phi.SSdU.a, base.phi.SSdU.b);
        phi_SSE = max(raw_SSE, phi_floor);
        phi_SSdU = max(raw_SSdU, phi_floor);
        phi_floored = (raw_SSE < phi_floor) || (raw_SSdU < phi_floor);
        phi_vintage = base.phi.vintage;
    else
        phi_SSE = 1;
        phi_SSdU = 1;
    end
    wall_phi_s = toc(t_phi);

    %% Controller
    t_build = tic;
    NMPC = build_nmpc(base, cfg, ...
        terminal_cost = opts.terminal_cost, ...
        set_setpoint = opts.set_setpoint, ...
        max_iter = opts.max_iter);
    wall_build_s = toc(t_build);

    %% Output skeleton
    out = struct();
    out.theta = theta(:).';
    out.cfg = cfg;
    out.tf = tf;
    out.N = N;
    out.T = base.T(1:N);
    out.extrapolated = opts.extrapolate;
    out.phi_SSE = phi_SSE;
    out.phi_SSdU = phi_SSdU;
    out.phi_floor = phi_floor;
    out.phi_floored = phi_floored;
    out.phi_vintage = phi_vintage;
    if opts.extrapolate
        out.phi_params_SSE = base.phi.SSE;
        out.phi_params_SSdU = base.phi.SSdU;
    else
        out.phi_params_SSE = [];
        out.phi_params_SSdU = [];
    end

    % Wall time of each stage. runtime_s below sums the solver time of the
    % control steps. These fields cover the rest of the call.
    out.wall_s = struct("phi", wall_phi_s, "build_nmpc", wall_build_s, ...
        "cases", 0, "total", 0);

    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;
    out.case = struct([]);

    n_cases = size(opts.x0, 1);

    %% Resume
    start_case = 1;
    resume_state = struct();
    if strlength(opts.checkpoint_path) > 0
        partial = load_partial_state(opts.checkpoint_path, opts.run_id);
        if partial.valid
            out = partial.out;
            start_case = partial.case_id;
            resume_state = partial.case_state;
            fprintf("Resuming from partial checkpoint (case %d): %s\n", start_case, opts.checkpoint_path);
        end
    end

    if start_case > n_cases
        out = aggregate_cases(out);
        return
    end

    %% Cases
    t_cases = tic;
    for case_id = start_case:n_cases
        if case_id == start_case && ~isempty(fieldnames(resume_state))
            case_resume = resume_state;
        else
            case_resume = struct();
        end

        if strlength(opts.checkpoint_path) > 0
            checkpoint_fn = @(case_state, i) save_case_checkpoint( ...
                opts.checkpoint_path, opts.run_id, out, case_id, case_state, ...
                base, opts.x0(case_id, :), i);
        else
            checkpoint_fn = [];
        end

        case_out = nmpc_run_case(base, NMPC, N, opts.x0(case_id, :), case_id, ...
            setpoint_fn = opts.setpoint_fn, ...
            checkpoint_fn = checkpoint_fn, ...
            resume = case_resume, ...
            verbosity = opts.verbosity, ...
            run_id = opts.run_id, ...
            log_path = opts.log_path);

        % The costs measured at the simulated fidelity survive the scaling.
        case_out.SSE_measured = case_out.SSE;
        case_out.SSdU_measured = case_out.SSdU;
        case_out.SSE = case_out.SSE_measured / phi_SSE;
        case_out.SSdU = case_out.SSdU_measured / phi_SSdU;
        case_out.cost_total = case_out.SSE + 1e4 * case_out.SSdU;

        if isempty(out.case)
            out.case = case_out;
        else
            out.case(case_id) = case_out;
        end

        out = aggregate_cases(out);

        if strlength(opts.checkpoint_path) > 0
            save_partial_state(opts.checkpoint_path, opts.run_id, out, case_id + 1, struct());
        end
    end

    out = aggregate_cases(out);
    out.wall_s.cases = toc(t_cases);

    % Flush any line that the logger had to buffer because the file was busy.
    if strlength(opts.log_path) > 0
        log_simulation_event(opts.log_path);
    end

    out.wall_s.total = toc(t_call);
end

function out = aggregate_cases(out)
%AGGREGATE_CASES Sum the per-case costs and runtimes into the run totals.
% Both the scaled and the measured costs are summed. A case summarised at a
% mid-run checkpoint carries only the measured value, because finalize_case
% runs before the scaling, so the measured total falls back to that field.
    out.SSE = 0;
    out.SSdU = 0;
    out.SSE_measured = 0;
    out.SSdU_measured = 0;
    out.runtime_s = 0;
    out.n_flag_not_one = 0;
    if ~isfield(out, "case") || isempty(out.case) || isempty(fieldnames(out.case))
        return
    end
    for k = 1:numel(out.case)
        c = out.case(k);
        if isfield(c, "SSE") && isfinite(c.SSE), out.SSE = out.SSE + c.SSE; end
        if isfield(c, "SSdU") && isfinite(c.SSdU), out.SSdU = out.SSdU + c.SSdU; end
        if isfield(c, "runtime_s") && isfinite(c.runtime_s), out.runtime_s = out.runtime_s + c.runtime_s; end
        if isfield(c, "n_flag_not_one"), out.n_flag_not_one = out.n_flag_not_one + c.n_flag_not_one; end

        if isfield(c, "SSE_measured")
            m_SSE = c.SSE_measured;
        elseif isfield(c, "SSE")
            m_SSE = c.SSE;
        else
            m_SSE = NaN;
        end
        if isfield(c, "SSdU_measured")
            m_SSdU = c.SSdU_measured;
        elseif isfield(c, "SSdU")
            m_SSdU = c.SSdU;
        else
            m_SSdU = NaN;
        end
        if isfinite(m_SSE), out.SSE_measured = out.SSE_measured + m_SSE; end
        if isfinite(m_SSdU), out.SSdU_measured = out.SSdU_measured + m_SSdU; end
    end
    out.J = out.SSE + 1e4 * out.SSdU;
end

function save_case_checkpoint(path, run_id, out, case_id, case_state, base, x0, i)
%SAVE_CASE_CHECKPOINT Store a mid-case checkpoint with costs summarised so far.
    partial_case = finalize_case(base, x0, case_id, base.T(1:size(case_state.Y,1)), ...
        base.noise(1:size(case_state.Y,1), :), case_state.Y, case_state.Y_meas, ...
        case_state.Ysp, case_state.U, case_state.RUNTIME, case_state.EXITFLAG, i);

    if ~isfield(out, "case") || isempty(out.case)
        out.case = partial_case;
    else
        out.case(case_id) = partial_case;
    end
    out = aggregate_cases(out);
    save_partial_state(path, run_id, out, case_id, case_state);
end

function save_partial_state(path, run_id, out, case_id, case_state)
%SAVE_PARTIAL_STATE Write the resume file for an interrupted run.
    partial = struct();
    partial.run_id = string(run_id);
    partial.case_id = case_id;
    partial.out = out;
    partial.case_state = case_state;
    partial.saved_at = char(datetime("now"));
    save(path, "partial");
end

function partial = load_partial_state(path, run_id)
%LOAD_PARTIAL_STATE Read a resume file, rejecting one written for another run.
    partial = struct("valid", false, "case_id", 1, "out", struct(), "case_state", struct());
    if ~isfile(path)
        return
    end
    S = load(path, "partial");
    if ~isfield(S, "partial")
        warning("Ignoring malformed partial file (missing struct): %s", path);
        return
    end
    p = S.partial;
    if ~isfield(p, "run_id") || string(p.run_id) ~= string(run_id)
        warning("Ignoring partial file written for a different run id: %s", path);
        return
    end
    if ~isfield(p, "case_id") || ~isfield(p, "out") || ~isfield(p, "case_state")
        warning("Ignoring malformed partial file (missing required fields): %s", path);
        return
    end
    partial.valid = true;
    partial.case_id = double(p.case_id);
    partial.out = p.out;
    partial.case_state = p.case_state;
end
