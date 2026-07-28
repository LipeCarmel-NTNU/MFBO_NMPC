function J = tuning_cost(var, f, plant, xss, uss, scenarios, opts)
    % TUNING_COST  Nonlinear closed-loop cost of an LQR tuning vector.
    %
    %   J = tuning_cost(var, f, plant, xss, uss, scenarios)
    %   J = tuning_cost(var, ..., y_weight = [10 1 1], reg_du = 1e4, ...)
    %
    % For the candidate tuning vector, the weight map builds (Q, R, ...), a gain
    % is synthesized with lqr_synth, and the nonlinear closed loop is simulated
    % from every initial condition in `scenarios` (cell array of x0). J is the
    % mean scenario performance plus regularization. Intended usage:
    %
    %   obj = @(var) tuning_cost(var, f, plant, xss, uss, scenarios, ...);
    %   [var_opt, J_opt] = lqr_tune(obj, var0);
    %
    % Name-value options
    %   weight_map : @(var, plant) weights struct, default @log10_weights
    %   Tf         : simulation horizon per scenario, default 20
    %   y_weight   : state-error weights in the performance index, default ones
    %   perf       : @(Y, T, U, xss) scalar; overrides the default weighted SSE
    %   reg_du     : coefficient on sum of squared input moves, default 0
    %   reg_var    : coefficient on ||var||^2, default 0
    %   u_min, u_max, ode_opt : passed to lqr_sim

    arguments
        var double
        f function_handle
        plant struct
        xss double
        uss double
        scenarios
        opts.weight_map function_handle = @log10_weights
        opts.Tf double = 20
        opts.y_weight double = ones(1, plant.nx)
        opts.perf = []
        opts.reg_du double = 0
        opts.reg_var double = 0
        opts.u_min double = -inf(plant.nu, 1)
        opts.u_max double =  inf(plant.nu, 1)
        opts.ode_opt = odeset()
    end

    if isempty(opts.perf)
        opts.perf = @(Y, T, U, xref) default_perf(Y, U, xref, opts.y_weight, opts.reg_du);
    end
    if ~iscell(scenarios)
        scenarios = {scenarios};
    end

    % Synthesize the gain; infeasible weight combinations get a large penalty
    try
        weights = opts.weight_map(var, plant);
        K = lqr_synth(plant, weights);
    catch
        J = 1e10;
        return
    end

    J = 0;
    for i = 1:numel(scenarios)
        [Y, T, U] = lqr_sim(f, plant, K, scenarios{i}, xss, uss, opts.Tf, ...
            u_min = opts.u_min, u_max = opts.u_max, ode_opt = opts.ode_opt);
        J = J + opts.perf(Y, T, U, xss);
    end
    J = J / numel(scenarios);

    % Keep tuning parameters from drifting to extreme magnitudes
    J = J + opts.reg_var * sum(var(:).^2);
end

function J = default_perf(Y, U, xref, y_weight, reg_du)
    % Weighted sum of squared state errors plus input-move regularization
    E = (Y - xref(:).') .* y_weight(:).';
    J = sum(E(:).^2) + reg_du * sum(sum(diff(U).^2));
end
