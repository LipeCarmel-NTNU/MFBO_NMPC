function [P, K, Ai, Bi, xss, uss, LQR_data] = terminal_P_xu_du( ...
        Vsp, Xsp, par, model, ode_opt, Ts, log10w, Q_eval, R1_eval, R2_eval)
%TERMINAL_P_XU_DU Terminal cost for the incremental LQR with stage cost
%x'Qx + u'R1u + du'R2du.
%
%   The steady state for the requested setpoint is linearised and converted to
%   incremental form, an LQR gain is designed for the fixed policy du = -K z,
%   and the terminal matrix P is the infinite-horizon cost of that policy
%   evaluated with the supplied weights.
%
%   z = [x - xss; u_prev - uss], so
%       x = Sx z
%       u - uss = (Su - K) z

    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);

    [A, B] = linearize(xss, uss, model);
    [Ai, Bi] = incremental(A, B, Ts);

    nx = size(A, 1);
    nu = size(B, 2);

    [K, ~, ~, ~] = build_LQR_full(log10w, Ai, Bi, nx, nu);

    % Closed loop under the fixed incremental policy du = -K z.
    Acl = Ai - Bi * K;

    Sx = [eye(nx), zeros(nx, nu)];
    Su = [zeros(nu, nx), eye(nu)];

    LQR_data = struct("Sx", Sx, "K", K, "Acl", Acl, "Su", Su);
    P = construct_P(LQR_data, Q_eval, R1_eval, R2_eval);
end
