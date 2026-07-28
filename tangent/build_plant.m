function plant = build_plant(A, B, Ts, mode)
    % BUILD_PLANT  Discretize (A, B) and package the design model for the suite.
    %
    %   plant = build_plant(A, B, Ts, 'positional')
    %   plant = build_plant(A, B, Ts, 'incremental')
    %
    % Positional:  x_{k+1} = Ad x_k + Bd u_k, design variable is u - uss.
    % Incremental: z = [x - xss; u_{k-1} - uss], design variable is du, with
    %
    %   z_{k+1} = Ai z_k + Bi du_k,
    %   Ai = [Ad, Bd; 0, I],  Bi = [Bd; I].
    %
    % Fields: mode, Ts, nx, nu, Ad, Bd, and for incremental additionally
    % Ai, Bi, Sx (z -> x-deviation), Su (z -> previous-input deviation).

    nx = size(A, 1);
    nu = size(B, 2);

    sysd = c2d(ss(A, B, eye(nx), zeros(nx, nu)), Ts, 'zoh');
    Ad = sysd.A;
    Bd = sysd.B;

    plant = struct('mode', mode, 'Ts', Ts, 'nx', nx, 'nu', nu, 'Ad', Ad, 'Bd', Bd);

    switch mode
        case 'positional'
            % Nothing further
        case 'incremental'
            plant.Ai = [Ad, Bd; zeros(nu, nx), eye(nu)];
            plant.Bi = [Bd; eye(nu)];
            plant.Sx = [eye(nx), zeros(nx, nu)];
            plant.Su = [zeros(nu, nx), eye(nu)];
        otherwise
            error('build_plant:mode', 'mode must be ''positional'' or ''incremental''.')
    end
end
