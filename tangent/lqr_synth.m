function [K, design] = lqr_synth(plant, weights)
    % LQR_SYNTH  Synthesize a discrete LQR gain for a positional or incremental plant.
    %
    %   [K, design] = lqr_synth(plant, weights)
    %
    % Positional (plant.mode = 'positional'):
    %   weights.Q (nx x nx), weights.R (nu x nu)
    %   Stage cost x'Qx + u'Ru in deviation variables; u = uss - K (x - xss).
    %
    % Incremental (plant.mode = 'incremental'):
    %   weights.Q (nx x nx), weights.R1 (nu x nu), weights.R2 (nu x nu)
    %   Stage cost x'Qx + u'R1u + du'R2du with u = u_{k-1} + du. In the
    %   augmented state z this expands to the cross-term form
    %     z'Qz z + 2 z'N du + du'R du,
    %   with Qz = blkdiag(Q, R1), N = [0; R1], R = R1 + R2, and du = -K z.
    %
    % design returns the matrices actually passed to dlqr.

    switch plant.mode
        case 'positional'
            design.Qz = weights.Q;
            design.R  = weights.R;
            K = dlqr(plant.Ad, plant.Bd, design.Qz, design.R);
            design.Acl = plant.Ad - plant.Bd * K;

        case 'incremental'
            design.Qz = blkdiag(weights.Q, weights.R1);
            design.R  = weights.R1 + weights.R2;
            design.N  = [zeros(plant.nx, plant.nu); weights.R1];
            K = dlqr(plant.Ai, plant.Bi, design.Qz, design.R, design.N);
            design.Acl = plant.Ai - plant.Bi * K;

        otherwise
            error('lqr_synth:mode', 'Unknown plant mode ''%s''.', plant.mode)
    end
end
