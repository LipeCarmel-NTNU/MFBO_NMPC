function [xss, uss, info] = trim_point(f, yss, x_guess, uss, opts)
    % TRIM_POINT  Find a steady state of xdot = f(x,u) subject to output targets.
    %
    %   [xss, uss, info] = trim_point(f, yss, x_guess, uss)
    %   [xss, uss, info] = trim_point(f, yss, x_guess, uss, h = @(x) x(2))
    %
    % NaN entries mark free quantities:
    %   yss(i) = NaN : output i is unconstrained
    %   uss(i) = NaN : input i is free, solved for (initial guess u_guess(i))
    %   uss(i) = val : input i is fixed at val
    %
    % With the default h(x) = x, constraining y(i) fixes state i directly, so
    % mixed fixed/free states are expressed through yss.
    %
    % The system [f(x,u); h(x)(cy) - yss(cy)] = 0, cy = ~isnan(yss), has
    % nx + #constrained_y equations in nx + #free_u unknowns. Levenberg-
    % Marquardt handles non-square systems in a least-squares sense; check
    % info.resnorm regardless.

    arguments
        f function_handle
        yss double
        x_guess double
        uss double
        opts.h function_handle = @(x) x(:)
        opts.u_guess double = zeros(size(uss))
        opts.options = []
    end

    h = opts.h;

    x_guess = x_guess(:);
    yss     = yss(:);
    uss     = uss(:);
    nx = numel(x_guess);

    cy     = ~isnan(yss);  % constrained outputs
    free_u =  isnan(uss);  % free inputs

    % Working input vector: fixed entries from uss, free entries from the guess
    u_work = uss;
    u_guess = opts.u_guess(:);
    u_work(free_u) = u_guess(free_u);

    if isempty(opts.options)
        opts.options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
            'Display', 'off');
    end

    % Solve for [x; u(free)]
    var0 = [x_guess; u_work(free_u)];
    [var, ~, exitflag] = fsolve(@(v) residual(v), var0, opts.options);

    xss = var(1:nx);
    uss = u_work;
    uss(free_u) = var(nx+1:end);

    res = residual(var);
    info.exitflag     = exitflag;
    info.resnorm      = norm(res);
    info.residual     = res;
    info.n_equations  = nx + nnz(cy);
    info.n_unknowns   = nx + nnz(free_u);
    if info.resnorm > 1e-6
        warning('trim_point:residual', ...
            'Steady-state residual norm %.3g exceeds 1e-6. Check the trim point.', info.resnorm)
    end

    function r = residual(v)
        x = v(1:nx);
        u = u_work;
        u(free_u) = v(nx+1:end);
        y = h(x);
        r = [f(x, u); y(cy) - yss(cy)];
    end
end
