function [K, Qi, Ri, Ai, Bi] = build_LQR(log10w, A, B, Ts)
    w = 10.^log10w;
    q = w(1:2);
    r = w(3:end);

    Q = diag([1 q]);
    R = diag(r);

    % LQR
    n = size(A,1);
    m = size(B,2);

    % Discretise
    sysc = ss(A,B,eye(n),zeros(n,m));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    % incremental form
    Ai = [Ad, Bd;
        zeros(m,n), eye(m)];
    Bi = [Bd;
        eye(m)];


    Qi = blkdiag(Q, zeros(m));  % Q only on the natural states
    Ri = R;                     % R acts on delta u explicitly

    try
        K = dlqr(Ai, Bi, Qi, Ri);
    catch
        warning('something went wrong. Probably infinite cost somewhere')
        keyboard
    end

    % The following yields an interesting structure, but cant be trusted
    % Ai = [A, B;
    %       zeros(m,n), zeros(m)];
    % Bi = [B;
    %       eye(m)];
    % Qi = blkdiag(Q, zeros(m));   % penalise x only
    % Ri = R*Ts^2;                 % The input is dudt and scalling matters
    % [K2,S,e] = lqrd(Ai,Bi,Qi,Ri,Ts);

end