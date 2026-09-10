function gen_data_timeline(ctx)
%GEN_DATA_TIMELINE Build one DOE-then-BO execution timeline per case.
%
%   Partial preprocessing, not a result: writes storage/timeline.mat with, per
%   case, the evaluations in execution order carrying the derived timing
%   columns the runtime diagnostics need:
%     totalH  wall time of the evaluation, hours
%     nmpcH   NMPC solve time, hours
%     gapMin  totalH - nmpcH, minutes
%     nBad    solves whose fmincon flag was 0 or negative
%
%   The figure and the numerical summary both read this, so the derivation
%   lives in one place and the two cannot disagree.

    P     = ra_require(ctx, "preprocessed");
    evals = P.evals;
    doe   = P.doe;

    evals.case = removecats(categorical(evals.case));
    doe.case   = removecats(categorical(doe.case));
    caseNames  = categories(evals.case);
    caseNames  = caseNames(ra_order_cases(caseNames));   % MF cases first, baseline last
    nCases     = numel(caseNames);
    caseLabels = arrayfun(@pretty_case, string(caseNames));

    A        = cell(nCases, 1);
    doeCount = zeros(nCases, 1);
    for k = 1:nCases
        Dk = sortrows(doe(doe.case == caseNames{k},     :), 'iter');
        Ek = sortrows(evals(evals.case == caseNames{k}, :), 'iter');
        doeCount(k) = height(Dk);

        T = [Dk; Ek];
        T.phase_label = [repmat("DOE", height(Dk), 1); repmat("BO", height(Ek), 1)];
        T.gidx    = (1:height(T))';
        T.totalH  = double(T.t_total) / 3600;
        T.nmpcH   = double(T.t_nmpc)  / 3600;
        T.gapMin  = (double(T.t_total) - double(T.t_nmpc)) / 60;
        T.nBad    = double(T.n_flag0) + double(T.n_flag_neg);
        A{k} = T;
    end

    outPath = fullfile(ctx.storageDir, 'timeline.mat');
    save(outPath, 'A', 'doeCount', 'caseNames', 'caseLabels', 'nCases');
    fprintf('Wrote %s\n', outPath);
end
