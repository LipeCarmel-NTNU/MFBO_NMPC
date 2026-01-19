function [Y, T, U] = LQR_simulation(system, Ts, Tf, y0, K, yss, uss, ode_opt)

    % Steps
    num_sim = ceil(Tf/Ts);

    % Dimensions
    n = length(yss);
    m = length(uss);

    % Preallocate
    Y = zeros(num_sim, n);
    U = zeros(num_sim, m);
    T = zeros(num_sim, 1);

    % Initialise
    Y(1,:) = y0(:).';
    T(1)   = 0;

    y_current = y0(:);
    t_current = 0;

    % Initialise previous input at setpoint (piecewise-constant hold)
    u_prev = uss(:);

    U(1,:) = u_prev(:).';

    for k = 2:num_sim
        % Deviations
        x_tilde      = y_current - yss(:);
        u_tilde_prev = u_prev    - uss(:);

        % Augmented deviation state
        z = [x_tilde; u_tilde_prev];

        % Incremental control in deviation coordinates
        du_tilde = -K * z;

        % Absolute input update
        uk = u_prev + du_tilde;

        % Constraints
        uk = max(uk, zeros(m,1));
        uk = min(uk, 0.4*ones(m,1));

        % ODE over the interval with piecewise-constant uk
        ode_current = @(t, y) system(t, y, uk);
        [t_span, y_span] = ode45(ode_current, [t_current, t_current + Ts], y_current, ode_opt);

        % Update state/time
        y_current = y_span(end,:).';
        t_current = t_span(end);

        % Store
        Y(k,:) = y_current.';
        U(k,:) = uk.';
        T(k)   = t_current;

        % Update previous input for next increment
        u_prev = uk;
    end
end
