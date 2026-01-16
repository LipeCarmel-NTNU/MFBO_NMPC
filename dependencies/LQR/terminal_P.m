function [P, K, Ai, Bi, xss, uss] = terminal_P(Vsp, Xsp, par, model, ode_opt, Ts, lqr_tuning, Q_eval, R_eval)
    % terminal_P  Infinite-horizon terminal matrix for a fixed incremental LQR policy.

    % Steady-state at setpoint
    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);

    % Linearise plant at (xss, uss)
    [A, B] = linearize(xss, uss, model);
    nu = size(B,2);

    % Incremental LQR gain and corresponding discrete-time lifted model
    [K, Qi, Ri, Ai, Bi] = build_LQR(lqr_tuning, A, B, Ts);

    Acl = Ai - Bi*K;
    Qz   = blkdiag(Q_eval, zeros(nu));
    Qbar = Qz + K' * R_eval * K;

    % Discrete-time infinite-horizon evaluation matrix
    P = dlyap(Acl', Qbar);
end
