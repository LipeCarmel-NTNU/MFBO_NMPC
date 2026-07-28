function [A, B, C] = linearize_fd(f, xss, uss, opts)
    % LINEARIZE_FD  Numeric linearization of xdot = f(x,u), y = h(x) at (xss, uss).
    %
    %   [A, B]    = linearize_fd(f, xss, uss)
    %   [A, B, C] = linearize_fd(f, xss, uss, h = @(x) x(1:2))
    %
    % Central finite differences with step eps^(1/3) * max(1, |v_i|). Numeric
    % differentiation is used instead of symbolic so f may contain interpolants,
    % lookups, or other non-symbolic constructs.
    %
    % Outputs
    %   A : df/dx at (xss, uss)   (nx x nx)
    %   B : df/du at (xss, uss)   (nx x nu)
    %   C : dh/dx at xss          (ny x nx), only if requested

    arguments
        f function_handle
        xss double
        uss double
        opts.h function_handle = @(x) x(:)
    end

    h = opts.h;

    xss = xss(:);
    uss = uss(:);

    A = jac(@(x) f(x, uss), xss);
    B = jac(@(u) f(xss, u), uss);
    if nargout > 2
        C = jac(h, xss);
    end
end

function J = jac(g, v0)
    n  = numel(v0);
    g0 = g(v0);
    m  = numel(g0);
    J  = zeros(m, n);
    for i = 1:n
        step = eps^(1/3) * max(1, abs(v0(i)));
        vp = v0; vp(i) = vp(i) + step;
        vm = v0; vm(i) = vm(i) - step;
        J(:, i) = (g(vp) - g(vm)) / (2*step);
    end
end
