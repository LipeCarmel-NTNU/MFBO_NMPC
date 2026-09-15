function gen_data_frontier(ctx)
%GEN_DATA_FRONTIER Split the evaluations per case and mark the Pareto sets.
%
%   Partial preprocessing, not a result: writes storage/frontier.mat with the
%   per-case BO tables (E), DOE tables (D), Pareto subsets (Tp), the boolean
%   masks (isPareto), the case names and labels, and the shared log axis
%   limits every objective-space figure draws to. Five figures and three
%   tables read this, so the masks are computed once.
%
%   Case membership comes from cases.txt via initial_preprocessing, so adding
%   a campaign is a data change, not a code change.

    P     = ra_require(ctx, "preprocessed");
    evals = P.evals;
    doe   = P.doe;

    evals.case = removecats(categorical(evals.case));
    doe.case   = removecats(categorical(doe.case));
    caseNames  = categories(evals.case);
    caseNames  = caseNames(ra_order_cases(caseNames));   % studied campaigns first, baseline last
    nCases     = numel(caseNames);
    if nCases > numel(ctx.caseMarkers)
        error('gen_data_frontier:tooManyCases', ...
            'Only %d case styles defined in ra_context; found %d cases.', ...
            numel(ctx.caseMarkers), nCases);
    end
    caseLabels = arrayfun(@pretty_case, string(caseNames));

    E = cell(nCases, 1); D = cell(nCases, 1);
    Tp = cell(nCases, 1); isPareto = cell(nCases, 1);

    for k = 1:nCases
        E{k} = sortrows(evals(evals.case == caseNames{k}, :), 'iter');
        D{k} = sortrows(doe(doe.case == caseNames{k},     :), 'iter');
        if isempty(E{k}) || isempty(D{k})
            error('gen_data_frontier:emptyCase', ...
                '%s has %d DOE and %d BO rows; both phases are required.', ...
                caseLabels(k), height(D{k}), height(E{k}));
        end

        isPareto{k} = compute_pareto_mask(double(E{k}.SSE), double(E{k}.SSdU));
        Tp{k}       = E{k}(isPareto{k}, :);
        [~, ord]    = sort(double(Tp{k}.SSE), 'descend');
        Tp{k}       = Tp{k}(ord, :);

        % Budgets differ per campaign on purpose: a run stopped on a wall-clock
        % budget has fewer BO rows than one stopped on an iteration count.
        fprintf('%-28s %3d DOE + %3d BO evaluations, %2d on its frontier\n', ...
            caseLabels(k), height(D{k}), height(E{k}), height(Tp{k}));
    end

    xLimAll = padded_log_limits(double(evals.SSdU), 1.25);
    yLimAll = padded_log_limits(double(evals.SSE),  1.25);
    % The colour axis spans the POOLED data, so the panels of one figure
    % share a scale and a marker means the same colour in each. A floor of 0.5
    % used to be forced here; with z running 0.833 to 1 that spent two thirds
    % of the ramp on an empty range and every marker came out one colour.
    zLo     = min(double(evals.z));
    zHi     = max(double(evals.z));
    if ~(zHi > zLo)
        % A campaign set can be all-z=1 (the single-fidelity arm on its own).
        % caxis needs an increasing pair, so open a nominal window below it.
        zLo = zHi - 0.01;
    end

    outPath = fullfile(ctx.storageDir, 'frontier.mat');
    save(outPath, 'E', 'D', 'Tp', 'isPareto', 'caseNames', 'caseLabels', ...
        'nCases', 'xLimAll', 'yLimAll', 'zLo', 'zHi');
    fprintf('Wrote %s\n', outPath);
end
