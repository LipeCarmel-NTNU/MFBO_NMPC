function case_out = nmpc_run_case(base, NMPC, N, x0, case_id, opts)
%NMPC_RUN_CASE Closed-loop simulation of one initial-condition case.
%
%   case_out = nmpc_run_case(base, NMPC, N, x0, case_id) runs N control steps
%   from x0, applying the measurement noise held in base.noise and integrating
%   the plant with ode45 between steps.
%
%   NMPC is a handle object, so the warm start it accumulates in latest_wopt
%   persists across the steps of this case and into any case run afterwards
%   with the same object.
%
%   Name-value options:
%     setpoint_fn   function handle called before each step as
%                   sp_state = setpoint_fn(NMPC, t_h, xk, case_id, sp_state,
%                   i, is_first), for runs whose setpoint or terminal cost
%                   changes over time. Called outside the iteration timer so
%                   that recomputing a terminal matrix is not counted as solve
%                   time. sp_state carries the hook's own state between steps
%                   and is stored in case_out.sp_state.
%     checkpoint_fn function handle called as checkpoint_fn(case_state, i)
%                   every ckpt_every_s seconds and on the last step, for long
%                   runs that must survive an interruption. The seconds it
%                   spends are returned in case_out.wall_ckpt_s and are not in
%                   case_out.runtime_s.
%     resume        state returned in an earlier case_state, to continue a
%                   partially completed case. It carries the solver state,
%                   which is written back into NMPC, so the resumed steps warm
%                   start exactly as they would have without the interruption.
%     verbosity     "full" prints the measurement, the action and a progress
%                   estimate; "control" prints the action only; "quiet"
%                   prints nothing.
%     run_id        identifier written to the log alongside a non-ideal
%                   optimiser flag.
%     log_path      SIMULATIONS_LOG.txt path. Empty disables flag logging.
%
%   case_out records the trajectories, the per-step runtime, and the fmincon
%   exit flag of every step in EXITFLAG.

    arguments
        base struct
        NMPC
        N (1,1) double
        x0 (1,:) double
        case_id (1,1) double
        opts.setpoint_fn = []
        opts.checkpoint_fn = []
        opts.resume struct = struct()
        opts.verbosity (1,1) string {mustBeMember(opts.verbosity, ["full" "control" "quiet"])} = "full"
        opts.run_id (1,1) string = ""
        opts.log_path (1,1) string = ""
    end

    noise = base.noise(1:N, :);
    T = base.T(1:N);

    %% Logging arrays
    i_start = 1;
    uk = zeros(1, base.nu);
    Y = zeros(N, base.nx);
    Y_meas = zeros(N, base.nx);
    Ysp = zeros(N, base.nx);
    U = zeros(N, base.nu);
    RUNTIME = zeros(N, 1);
    EXITFLAG = nan(N, 1);
    xk = x0;

    % State carried by a setpoint_fn across steps, opaque to this function.
    sp_state = struct();

    if ~isempty(fieldnames(opts.resume))
        r = opts.resume;
        i_start = max(1, min(N, double(r.i_next)));
        if isfield(r, "uk"), uk = double(r.uk); end
        if isfield(r, "xk"), xk = double(r.xk); end
        if isfield(r, "Y"), Y = double(r.Y); end
        if isfield(r, "Y_meas"), Y_meas = double(r.Y_meas); end
        if isfield(r, "Ysp"), Ysp = double(r.Ysp); end
        if isfield(r, "U"), U = double(r.U); end
        if isfield(r, "RUNTIME"), RUNTIME = double(r.RUNTIME); end
        if isfield(r, "EXITFLAG"), EXITFLAG = double(r.EXITFLAG); end
        if isfield(r, "sp_state"), sp_state = r.sp_state; end
        % latest_wopt and latest_flag are the whole per-step solver state, and
        % solve() branches on both: a negative flag selects a cold guess and
        % the infeasibility fallback. NMPC is a handle object, so writing them
        % back here resumes the case from the state it had reached. Without
        % this the step after the interruption solves cold, which changes the
        % trajectory and inflates its runtime.
        if isfield(r, "wopt"), NMPC.latest_wopt = r.wopt; end
        if isfield(r, "flag"), NMPC.latest_flag = r.flag; end
    end

    has_setpoint_fn = ~isempty(opts.setpoint_fn);
    has_checkpoint_fn = ~isempty(opts.checkpoint_fn);
    do_log = strlength(opts.log_path) > 0;

    % Cap on the flag lines that one case may write. See the logging block in
    % the loop below.
    max_flag_log = 20;
    n_flagged = 0;

    % Checkpoints go on the wall clock and not on simulated time, so an
    % interruption costs at most this much work whatever the fidelity. The
    % whole hour of simulated time that this used to fire on never came around
    % at all in a run shorter than an hour.
    ckpt_every_s = 600;
    ckpt_timer = tic;
    wall_ckpt_s = 0;

    %% Control loop
    for i = i_start:N
        t_now = T(i);

        % Setpoint and terminal-cost updates happen before the timer starts.
        if has_setpoint_fn
            sp_state = opts.setpoint_fn(NMPC, t_now, xk, case_id, sp_state, i, i == i_start);
        end

        iter_timer = tic;

        Y(i, :) = xk;
        Ysp(i, :) = NMPC.x_sp(1:base.nx);

        % Measurement is the state plus Gaussian noise, clipped to non-negative.
        yk_meas = xk + noise(i, :);
        yk_meas(yk_meas < 0) = 0;
        Y_meas(i, :) = yk_meas;

        uk = NMPC.solve(yk_meas(:)', uk(:)');
        U(i, :) = uk;

        % NMPC.solve leaves the fmincon exit flag of this step in latest_flag.
        %
        % The log is capped per case. A flag other than 1 is meant to be the
        % exception, but a tight iteration budget makes every step report flag 0,
        % and one open-write-close for each control step then dominates the run
        % time on a synced folder. The cap keeps the first max_flag_log events
        % and writes one summary line at the end of the case. The full flag
        % vector is stored in the .mat either way, so nothing is lost.
        EXITFLAG(i) = NMPC.latest_flag;
        if do_log && ~(EXITFLAG(i) == 1)
            n_flagged = n_flagged + 1;
            if n_flagged <= max_flag_log
                log_simulation_event(opts.log_path, sprintf( ...
                    "%s | id=%s | case=%d | step=%d/%d | t=%.4f h | flag=%g", ...
                    char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
                    opts.run_id, case_id, i, N, t_now, EXITFLAG(i)));
            end
        end

        if opts.verbosity == "full"
            fprintf('\nFlag: \n')
            disp(NMPC.latest_flag)
            fprintf('\nMeasurement: \n')
            disp(yk_meas)
            fprintf('\nControl action: \n')
            disp(uk)
        elseif opts.verbosity == "control"
            disp(uk)
        end

        % One-step plant integration. The last sample needs no propagation.
        if N > 1 && i < N
            [~, y] = ode45(@(t, x) base.plant(x, uk), base.tspan, xk, base.ode_opt);
            xk = y(end, :);
        end

        RUNTIME(i) = toc(iter_timer);

        if opts.verbosity == "full"
            elapsed_min = sum(RUNTIME(1:i)) / 60;
            progress = i / max(N, 1);
            total_min = elapsed_min / max(progress, eps);
            fprintf('Case %d: %.1f %% | elapsed %.2f min | total %.2f min | left %.2f min\n', ...
                case_id, 100*progress, elapsed_min, total_min, total_min - elapsed_min);
        end

        % The test and the write sit after RUNTIME(i) is recorded, so neither
        % enters the per-step time that case_out.runtime_s sums.
        if has_checkpoint_fn && (i == N || toc(ckpt_timer) >= ckpt_every_s)
            case_state = struct();
            case_state.i_next = min(i + 1, N + 1);
            case_state.uk = uk;
            case_state.xk = xk;
            case_state.Y = Y;
            case_state.Y_meas = Y_meas;
            case_state.Ysp = Ysp;
            case_state.U = U;
            case_state.RUNTIME = RUNTIME;
            case_state.EXITFLAG = EXITFLAG;
            case_state.sp_state = sp_state;
            case_state.wopt = NMPC.latest_wopt;
            case_state.flag = NMPC.latest_flag;
            t_ck = tic;
            opts.checkpoint_fn(case_state, i);
            wall_ckpt_s = wall_ckpt_s + toc(t_ck);
            ckpt_timer = tic;
        end
    end

    % One summary line when the cap suppressed flag events, so the log states
    % how many steps it did not list.
    if do_log && n_flagged > max_flag_log
        log_simulation_event(opts.log_path, sprintf( ...
            "%s | id=%s | case=%d | %d of %d steps had a flag other than 1, " + ...
            "%d line(s) suppressed by the cap of %d", ...
            char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
            opts.run_id, case_id, n_flagged, N, n_flagged - max_flag_log, max_flag_log));
    end

    case_out = finalize_case(base, x0, case_id, T, noise, Y, Y_meas, Ysp, U, RUNTIME, EXITFLAG, N);

    % Wall time this case spent on checkpoints, for the caller to subtract.
    case_out.wall_ckpt_s = wall_ckpt_s;

    % Whatever the setpoint hook accumulated, for callers that log a schedule.
    case_out.sp_state = sp_state;
end
