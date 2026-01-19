function P = construct_P(LQR_data, Q, R1, R2)

    Sx  = LQR_data.Sx;
    Su  = LQR_data.Su;
    K   = LQR_data.K;
    Acl = LQR_data.Acl;

    Qbar = (Sx.'*Q*Sx) ...
         + (Su - K).'*R1*(Su - K) ...
         + K.'*R2*K;

    % Infinite-horizon evaluation matrix
    P = dlyap(Acl', Qbar);
end