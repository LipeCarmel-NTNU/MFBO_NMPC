function [K, Qz, R, N] = build_LQR_full(log10w, Ai, Bi, nx, nu)
    % Stage cost in (z, du):
    %   x'Qx + u'R1u + du'R2du, with u = u_{k-1} + du
    %
    % Expanded form:
    %   z'Qz z + 2 z' N du + du' R du
    %
    % with Qz = blkdiag(Q, R1), N = [0; R1], R = R1 + R2.

    w = 10.^log10w(:);

    nq = nx - 1; % Q(1,1) = 1

    q  = [1; w(1:nq)];
    r1 = w(nq + (1:nu));
    r2 = w(nq + nu + (1:nu));

    Q  = diag(q);
    R1 = diag(r1);
    R2 = diag(r2);

    Qz = blkdiag(Q, R1);
    R  = R1 + R2;
    N  = [zeros(nx,nu);
        R1];

    try
        K = dlqr(Ai, Bi, Qz, R, N);
    catch ME
        warning('Code is about to crash. See:')
        warning(ME.message)
        keyboard
        throw(ME)
    end
end
