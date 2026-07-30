function fn = make_schedule_setpoint_fn(base, schedule, is_benchmark, lqr_tuning)
%MAKE_SCHEDULE_SETPOINT_FN Setpoint hook for a piecewise-constant Xsp schedule.
%
%   fn = make_schedule_setpoint_fn(base, schedule, is_benchmark, lqr_tuning)
%   returns a handle for the setpoint_fn option of nmpc_run_case.
%
%   schedule.segment_end_h and schedule.Xsp_values define the piecewise
%   constant setpoint for the second state; schedule.Vsp is the first-state
%   setpoint. On the first step, and whenever the segment changes, the steady
%   state is recomputed and the terminal matrix is rebuilt: from find_ss with
%   P = 0 for the benchmark controller, and from terminal_P_xu_du otherwise.
%   The terminal reference x_term/u_term follows the setpoint on every step.
%
%   Case 2 additionally overwrites the third setpoint with
%   Ssp = min(3, 2*(Xsp - X)), which may be negative by design.

    fn = @(NMPC, t_h, xk, case_id, sp_state, i, is_first) apply( ...
        NMPC, t_h, xk, case_id, sp_state, i, is_first, base, schedule, is_benchmark, lqr_tuning);
end

function sp_state = apply(NMPC, t_h, xk, case_id, sp_state, i, is_first, base, schedule, is_benchmark, lqr_tuning)

    if ~isfield(sp_state, "currentXsp2")
        sp_state.currentXsp2 = NaN;
        sp_state.Xsp2_trace = zeros(base.N, 1);
    end

    xsp2_target = scheduled_xsp2(schedule, t_h);

    if is_first || abs(xsp2_target - sp_state.currentXsp2) > 1e-12
        if is_benchmark
            [xss, uss] = find_ss(schedule.Vsp, xsp2_target, base.par, base.model, base.ode_opt);
            P = zeros(base.nx + base.nu);
        else
            [P, ~, ~, ~, xss, uss, ~] = terminal_P_xu_du( ...
                schedule.Vsp, xsp2_target, base.par, base.model, base.ode_opt, base.dt, ...
                lqr_tuning, NMPC.Q, NMPC.R_u, NMPC.R_du);
        end
        NMPC.x_sp = xss;
        NMPC.u_sp = uss;
        NMPC.P = P;
        sp_state.currentXsp2 = xsp2_target;
    end

    % Case 2 tracks a sugar setpoint that follows the biomass gap.
    if case_id == 2 && numel(NMPC.x_sp) >= 3
        NMPC.x_sp(3) = min(3, 2 * (sp_state.currentXsp2 - xk(2)));
    end

    % P is augmented, so it needs the reference point of its own coordinates.
    % The terminal reference is the tracking setpoint, including the sugar
    % override of case 2.
    NMPC.x_term = NMPC.x_sp;
    NMPC.u_term = NMPC.u_sp;

    idx = min(max(i, 1), numel(sp_state.Xsp2_trace));
    sp_state.Xsp2_trace(idx) = sp_state.currentXsp2;
end

function xsp2 = scheduled_xsp2(schedule, t_h)
    if t_h < schedule.segment_end_h(1)
        xsp2 = schedule.Xsp_values(1);
    elseif t_h < schedule.segment_end_h(2)
        xsp2 = schedule.Xsp_values(2);
    else
        xsp2 = schedule.Xsp_values(3);
    end
end
