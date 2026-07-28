function [Y, T, U] = lqr_sim(f, plant, K, x0, xss, uss, Tf, opts)
    % LQR_SIM  Closed-loop simulation of the nonlinear plant under the LQR law.
    %
    %   [Y, T, U] = lqr_sim(f, plant, K, x0, xss, uss, Tf)
    %   [Y, T, U] = lqr_sim(..., u_min = umin, u_max = umax, ode_opt = odeset(...))
    %
    % The input is held piecewise constant over each sample interval plant.Ts and
    % the nonlinear ODE xdot = f(x,u) is integrated with ode45.
    %
    % Control law
    %   positional : u_k = uss - K (x_k - xss)
    %   incremental: du_k = -K [x_k - xss; u_{k-1} - uss],  u_k = u_{k-1} + du_k
    %
    % Saturation (if given) is applied to the absolute input. Note that for the
    % incremental law the previous input stored for the next step is the
    % saturated one, which provides anti-windup behavior.

    arguments
        f function_handle
        plant struct
        K double
        x0 double
        xss double
        uss double
        Tf double
        opts.u_min double = -inf(plant.nu, 1)
        opts.u_max double =  inf(plant.nu, 1)
        opts.ode_opt = odeset()
        opts.u_init double = uss(:) % previous input at k = 0 (incremental only)
    end

    u_min   = opts.u_min(:);
    u_max   = opts.u_max(:);
    ode_opt = opts.ode_opt;
    u_prev  = opts.u_init(:);

    Ts  = plant.Ts;
    num = ceil(Tf/Ts) + 1;

    xss = xss(:);
    uss = uss(:);

    Y = zeros(num, plant.nx);
    U = zeros(num, plant.nu);
    T = zeros(num, 1);

    x = x0(:);
    t = 0;
    Y(1,:) = x.';
    U(1,:) = u_prev.';

    for k = 2:num
        switch plant.mode
            case 'positional'
                u = uss - K * (x - xss);
            case 'incremental'
                z  = [x - xss; u_prev - uss];
                u  = u_prev - K * z;
            otherwise
                error('lqr_sim:mode', 'Unknown plant mode ''%s''.', plant.mode)
        end

        % Input saturation
        u = min(max(u, u_min), u_max);

        [ts, xs] = ode45(@(tt, xx) f(xx, u), [t, t + Ts], x, ode_opt);
        x = xs(end, :).';
        t = ts(end);

        Y(k,:) = x.';
        U(k,:) = u.';
        T(k)   = t;

        u_prev = u;
    end
end
