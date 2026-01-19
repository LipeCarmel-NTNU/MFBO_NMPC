function [xss, uss] = find_ss(V, X, par, model, ode_opt)
    % Starvation:
    % Solve for low flow rate steady-state
    Fin = X * par.Y_XSinv * par.kd / par.Sin; % approximate solution

    var = [Fin 0]; % solve ss for Fin and S
    options = optimoptions('fsolve', 'Algorithm','levenberg-marquardt');
    var_opt = fsolve(@(var) model([V, X, var(2)], [var(1), 0, var(1)]), var, options);
    Fin = var_opt(1);
    S = var_opt(2);

    uss = [Fin 0 Fin];
    u = uss;
    TSPAN = [0 10];
    [t, x] = ode15s(@(t,x) model(x, u), TSPAN, [V, X, S], ode_opt); % check
    xss = x(end, :);
    if abs( X - xss(2) ) > 1e-5
        warning('Check steady-state')
    end
end