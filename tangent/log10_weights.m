function weights = log10_weights(var, plant)
    % LOG10_WEIGHTS  Map a log10 tuning vector to diagonal LQR weight matrices.
    %
    %   weights = log10_weights(var, plant)
    %
    % Q(1,1) is fixed to 1 to remove the scale invariance of the LQR cost, so
    % the free parameters are:
    %
    %   positional : var = [q(2:nx), r(1:nu)],           length nx-1 + nu
    %                Q = diag([1, 10.^q]), R = diag(10.^r)
    %   incremental: var = [q(2:nx), r1(1:nu), r2(1:nu)], length nx-1 + 2nu
    %                Q = diag([1, 10.^q]), R1 = diag(10.^r1), R2 = diag(10.^r2)

    var = var(:);
    nx = plant.nx;
    nu = plant.nu;
    nq = nx - 1;

    weights.Q = diag([1; 10.^var(1:nq)]);

    switch plant.mode
        case 'positional'
            assert(numel(var) == nq + nu, 'log10_weights: expected %d parameters.', nq + nu)
            weights.R = diag(10.^var(nq + (1:nu)));
        case 'incremental'
            assert(numel(var) == nq + 2*nu, 'log10_weights: expected %d parameters.', nq + 2*nu)
            weights.R1 = diag(10.^var(nq + (1:nu)));
            weights.R2 = diag(10.^var(nq + nu + (1:nu)));
        otherwise
            error('log10_weights:mode', 'Unknown plant mode ''%s''.', plant.mode)
    end
end
