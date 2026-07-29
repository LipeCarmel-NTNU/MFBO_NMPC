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
%   Extrapolation
%     opts.extrapolate divides the measured partial costs by the fitted cost
%     fraction frac(f) = I_f(a, b) to estimate full-horizon totals, and
%     requires base.beta. It is off during initialisation, where out.SSE and
%     out.SSdU are the costs actually measured at the simulated fidelity, which
%     is what the surrogate is later fitted to.
%
%     The measured costs survive the scaling as out.SSE_measured and
%     out.SSdU_measured, and the divisors and the surrogate vintage are stored
%     alongside them, so the full-horizon estimate can be recomputed under any
%     later fit without rerunning the simulation.
%
%   Name-value options:
%     horizon        "fidelity" (default) or "full"
%     extrapolate    estimate full-horizon costs from base.beta (default false)
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

    %% Fidelity surrogate fractions
    % frac is the share of the full-horizon cost accumulated by fidelity f.
    % I_f(a, b) is already confined to [0, 1] and equals 1 at f = 1, so the
    % floor of 0.01 is the only guard needed. It caps the scaling at 100x so
    % that a tiny fraction cannot turn a small partial cost into an enormous
    % estimate; whether it was active is recorded per evaluation.
    frac_floor = 0.01;
    beta_vintage = NaN;
    frac_floored = false;

    if opts.extrapolate
        if isempty(base.beta)
            error("opts.extrapolate is true but base.beta is empty. Build nmpc_base with beta_coeffs_path set, or reload it with load_beta_coeffs.");
        end
        raw_SSE = beta_eval(cfg.f, base.beta.SSE.a, base.beta.SSE.b);
        raw_SSdU = beta_eval(cfg.f, base.beta.SSdU.a, base.beta.SSdU.b);
        frac_SSE = max(raw_SSE, frac_floor);
        frac_SSdU = max(raw_SSdU, frac_floor);
        frac_floored = (raw_SSE < frac_floor) || (raw_SSdU < frac_floor);
        beta_vintage = base.beta.vintage;
    else
        frac_SSE = 1;
        frac_SSdU = 1;
    end

    %% Controller
    NMPC = build_nmpc(base, cfg, ...
        terminal_cost = opts.terminal_cost, ...
        set_setpoint = opts.set_setpoint, ...
        max_iter = opts.max_iter);

    %% Output skeleton
    out = struct();
    out.theta = theta(:).';
    out.cfg = cfg;
    out.tf = tf;
    out.N = N;
    out.T = base.T(1:N);
    out.extrapolated = opts.extrapolate;
    out.frac_SSE = frac_SSE;
    out.frac_SSdU = frac_SSdU;
    out.frac_floor = frac_floor;
    out.frac_floored = frac_floored;
    out.beta_vintage = beta_vintage;
    if opts.extrapolate
        out.beta_SSE = base.beta.SSE;
        out.beta_SSdU = base.beta.SSdU;
    else
        out.beta_SSE = [];
        out.beta_SSdU = [];
    end
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

        % Costs measured at the simulated fidelity, kept regardless of scaling.
        case_out.SSE_measured = case_out.SSE;
        case_out.SSdU_measured = case_out.SSdU;
        case_out.SSE = case_out.SSE_measured / frac_SSE;
        case_out.SSdU = case_out.SSdU_measured / frac_SSdU;
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

    % Flush anything the logger had to buffer because the file was busy.
    if strlength(opts.log_path) > 0
        log_simulation_event(opts.log_path);
    end
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
