function [var_opt, J_opt] = lqr_tune(objective_func, var0, opts)
    % LQR_TUNE  Minimize a tuning objective over the LQR tuning vector.
    %
    %   [var_opt, J_opt] = lqr_tune(objective_func, var0)
    %   [var_opt, J_opt] = lqr_tune(objective_func, var0, solver = 'fminsearch')
    %
    % objective_func : @(var) scalar. Build one with tuning_cost, e.g.
    %
    %   obj = @(var) tuning_cost(var, f, plant, xss, uss, scenarios, ...);
    %   [var_opt, J_opt] = lqr_tune(obj, zeros(1, plant.nx-1 + 2*plant.nu));
    %
    % Recover the result with:
    %
    %   K = lqr_synth(plant, log10_weights(var_opt, plant));

    arguments
        objective_func function_handle
        var0 double
        opts.solver {mustBeMember(opts.solver, {'fminunc', 'fminsearch'})} = 'fminunc'
        opts.solver_options = []
    end

    switch opts.solver
        case 'fminunc'
            if isempty(opts.solver_options)
                opts.solver_options = optimoptions('fminunc', ...
                    'Display', 'iter-detailed', 'MaxFunctionEvaluations', 2000);
            end
            [var_opt, J_opt] = fminunc(objective_func, var0, opts.solver_options);
        case 'fminsearch'
            if isempty(opts.solver_options)
                opts.solver_options = optimset('Display', 'iter', 'MaxFunEvals', 2000);
            end
            [var_opt, J_opt] = fminsearch(objective_func, var0, opts.solver_options);
    end
end
