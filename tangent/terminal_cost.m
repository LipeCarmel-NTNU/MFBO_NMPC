function [P, info] = terminal_cost(plant, K, eval_weights, opts)
    % TERMINAL_COST  Infinite-horizon cost-to-go matrix under a fixed LQR policy.
    %
    %   [P, info] = terminal_cost(plant, K, eval_weights)
    %   [P, info] = terminal_cost(plant, K, eval_weights, validate = true)
    %
    % Computes P such that (for the linear closed loop) the infinite-horizon sum
    % of the EVALUATION stage cost under the fixed feedback K equals z'Pz. The
    % evaluation weights are in general different from the weights used to
    % synthesize K (e.g. the NMPC stage-cost weights), in which case z'Pz is an
    % upper bound on the optimal cost-to-go for that stage cost, and hence a
    % valid NMPC terminal cost.
    %
    % Positional:  eval_weights.Q, eval_weights.R
    %   Qbar = Q + K'RK,                        Acl = Ad - Bd K,  z = x - xss
    % Incremental: eval_weights.Q, eval_weights.R1, eval_weights.R2
    %   with x = Sx z, u-deviation = (Su - K) z, du = -K z:
    %   Qbar = Sx'Q Sx + (Su-K)'R1(Su-K) + K'R2 K,  Acl = Ai - Bi K,
    %   z = [x - xss; u_{k-1} - uss]
    %
    % P solves the discrete Lyapunov equation P = Acl' P Acl + Qbar.
    %
    % Name-value options
    %   validate : if true (default), compare z0'Pz0 against a truncated sum of
    %              the linear closed-loop cost from a random z0; results in
    %              info.J_sum, info.J_P, info.rel_err
    %   N        : number of steps in the validation sum, default 5000

    arguments
        plant struct
        K double
        eval_weights struct
        opts.validate logical = true
        opts.N double = 5000
    end

    do_validate = opts.validate;
    N = opts.N;

    switch plant.mode
        case 'positional'
            Acl  = plant.Ad - plant.Bd * K;
            Qbar = eval_weights.Q + K.' * eval_weights.R * K;
        case 'incremental'
            Sx = plant.Sx;
            Su = plant.Su;
            Acl  = plant.Ai - plant.Bi * K;
            Qbar = Sx.' * eval_weights.Q * Sx ...
                 + (Su - K).' * eval_weights.R1 * (Su - K) ...
                 + K.' * eval_weights.R2 * K;
        otherwise
            error('terminal_cost:mode', 'Unknown plant mode ''%s''.', plant.mode)
    end

    if max(abs(eig(Acl))) >= 1
        error('terminal_cost:unstable', ...
            'Closed loop is not Schur stable (spectral radius %.4g); P does not exist.', ...
            max(abs(eig(Acl))))
    end

    P = dlyap(Acl.', Qbar);

    info.Acl  = Acl;
    info.Qbar = Qbar;

    if do_validate
        rng_state = rng; rng(0); % reproducible check, restore afterwards
        z0 = randn(size(Acl, 1), 1);
        rng(rng_state);

        J_sum = 0;
        z = z0;
        for k = 1:N
            J_sum = J_sum + z.' * Qbar * z;
            z = Acl * z;
        end
        info.J_sum   = J_sum;
        info.J_P     = z0.' * P * z0;
        info.rel_err = abs(info.J_sum - info.J_P) / max(1, abs(info.J_P));
        if info.rel_err > 1e-6
            warning('terminal_cost:validate', ...
                'Truncated-sum check rel. error %.3g; increase N or check conditioning.', ...
                info.rel_err)
        end
    end
end
