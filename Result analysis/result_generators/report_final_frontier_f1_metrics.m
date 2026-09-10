function report_final_frontier_f1_metrics(E, caseNames, projectRoot, numericalFolder)
%REPORT_FINAL_FRONTIER_F1_METRICS Final Pareto f=1 table with noisy/noiseless metrics.
% Guarded: skips with a warning when the z=1 re-evaluation MAT folders are absent.
rootFolder = fullfile(projectRoot, 'results');
refinedRoot = fullfile(rootFolder, 'final_fidelity_same_noise');
if ~isfolder(refinedRoot)
    warning('main_pareto:noF1Metrics', ...
        'Final f=1 metrics section skipped: %s not found.', refinedRoot);
    return
end

nCases = numel(E);
Tall = table();
for k = 1:nCases
    Tk = E{k};
    Tk.run_key = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.ts = ts_string(Tk.timestamp);
    Tall = [Tall; Tk]; %#ok<AGROW>
end
if isempty(Tall)
    fprintf('Final Pareto f=1 metrics table: no BO rows available.\n');
    return
end

isFront = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(isFront, :);
if isempty(Tf)
    fprintf('Final Pareto f=1 metrics table: no Pareto rows available.\n');
    return
end

n = height(Tf);
noisy_Jtrack = nan(n,1); noisy_JTV = nan(n,1);
noiseless_Jtrack = nan(n,1); noiseless_JTV = nan(n,1);

noisy_settle_x1_h_c1 = nan(n,1); noisy_settle_x2_h_c1 = nan(n,1); noisy_settle_x3_h_c1 = nan(n,1);
noisy_settle_x1_h_c2 = nan(n,1); noisy_settle_x2_h_c2 = nan(n,1); noisy_settle_x3_h_c2 = nan(n,1);
noisy_settle_case_h_c1 = nan(n,1); noisy_settle_case_h_c2 = nan(n,1);
noisy_IAE_x1_c1 = nan(n,1); noisy_IAE_x2_c1 = nan(n,1); noisy_IAE_x3_c1 = nan(n,1);
noisy_IAE_x1_c2 = nan(n,1); noisy_IAE_x2_c2 = nan(n,1); noisy_IAE_x3_c2 = nan(n,1);
noisy_IAE_case_c1 = nan(n,1); noisy_IAE_case_c2 = nan(n,1);

noiseless_settle_x1_h_c1 = nan(n,1); noiseless_settle_x2_h_c1 = nan(n,1); noiseless_settle_x3_h_c1 = nan(n,1);
noiseless_settle_x1_h_c2 = nan(n,1); noiseless_settle_x2_h_c2 = nan(n,1); noiseless_settle_x3_h_c2 = nan(n,1);
noiseless_settle_case_h_c1 = nan(n,1); noiseless_settle_case_h_c2 = nan(n,1);
noiseless_IAE_x1_c1 = nan(n,1); noiseless_IAE_x2_c1 = nan(n,1); noiseless_IAE_x3_c1 = nan(n,1);
noiseless_IAE_x1_c2 = nan(n,1); noiseless_IAE_x2_c2 = nan(n,1); noiseless_IAE_x3_c2 = nan(n,1);
noiseless_IAE_case_c1 = nan(n,1); noiseless_IAE_case_c2 = nan(n,1);

nMissingNoisy = 0;
for i = 1:n
    run_key = string(Tf.run_key(i));
    ts = string(Tf.ts(i));
    noisyOutPath = fullfile(refinedRoot, run_key + "_full_f1_same_noise", "out_full_" + ts + ".mat");
    noiselessOutPath = fullfile(rootFolder, 'test_run', run_key + "_full_f1_no_noise", "out_full_" + ts + ".mat");
    if ~isfile(noisyOutPath)
        nMissingNoisy = nMissingNoisy + 1;
        warning('Missing z=1 refined MAT for Pareto controller: %s | %s (%s)', run_key, ts, noisyOutPath);
    end
    if ~isfile(noiselessOutPath)
        warning('Missing no-noise MAT for Pareto controller: %s | %s (%s)', run_key, ts, noiselessOutPath);
    end

    noisyM = compute_out_metrics_by_path(noisyOutPath, 0.05);
    noiselessM = compute_out_metrics_by_path(noiselessOutPath, 0.05);

    noisy_Jtrack(i) = noisyM.Jtrack;
    noisy_JTV(i) = noisyM.JTV;
    noiseless_Jtrack(i) = noiselessM.Jtrack;
    noiseless_JTV(i) = noiselessM.JTV;

    noisy_settle_x1_h_c1(i) = noisyM.settle_h(1,1); noisy_settle_x2_h_c1(i) = noisyM.settle_h(1,2); noisy_settle_x3_h_c1(i) = noisyM.settle_h(1,3);
    noisy_settle_x1_h_c2(i) = noisyM.settle_h(2,1); noisy_settle_x2_h_c2(i) = noisyM.settle_h(2,2); noisy_settle_x3_h_c2(i) = noisyM.settle_h(2,3);
    noisy_settle_case_h_c1(i) = noisyM.settle_case_h(1); noisy_settle_case_h_c2(i) = noisyM.settle_case_h(2);
    noisy_IAE_x1_c1(i) = noisyM.IAE(1,1); noisy_IAE_x2_c1(i) = noisyM.IAE(1,2); noisy_IAE_x3_c1(i) = noisyM.IAE(1,3);
    noisy_IAE_x1_c2(i) = noisyM.IAE(2,1); noisy_IAE_x2_c2(i) = noisyM.IAE(2,2); noisy_IAE_x3_c2(i) = noisyM.IAE(2,3);
    noisy_IAE_case_c1(i) = noisyM.IAE_case(1); noisy_IAE_case_c2(i) = noisyM.IAE_case(2);

    noiseless_settle_x1_h_c1(i) = noiselessM.settle_h(1,1); noiseless_settle_x2_h_c1(i) = noiselessM.settle_h(1,2); noiseless_settle_x3_h_c1(i) = noiselessM.settle_h(1,3);
    noiseless_settle_x1_h_c2(i) = noiselessM.settle_h(2,1); noiseless_settle_x2_h_c2(i) = noiselessM.settle_h(2,2); noiseless_settle_x3_h_c2(i) = noiselessM.settle_h(2,3);
    noiseless_settle_case_h_c1(i) = noiselessM.settle_case_h(1); noiseless_settle_case_h_c2(i) = noiselessM.settle_case_h(2);
    noiseless_IAE_x1_c1(i) = noiselessM.IAE(1,1); noiseless_IAE_x2_c1(i) = noiselessM.IAE(1,2); noiseless_IAE_x3_c1(i) = noiselessM.IAE(1,3);
    noiseless_IAE_x1_c2(i) = noiselessM.IAE(2,1); noiseless_IAE_x2_c2(i) = noiselessM.IAE(2,2); noiseless_IAE_x3_c2(i) = noiselessM.IAE(2,3);
    noiseless_IAE_case_c1(i) = noiselessM.IAE_case(1); noiseless_IAE_case_c2(i) = noiselessM.IAE_case(2);
end

if nMissingNoisy == n
    warning('main_pareto:noF1Data', ...
        'Final f=1 metrics section skipped: no z=1 re-evaluation MATs found under %s.', refinedRoot);
    return
end

Treport = table( ...
    string(Tf.run_key), string(Tf.ts), double(Tf.Np), double(Tf.Nc), ...
    double(Tf.Q1), double(Tf.Q2), double(Tf.Q3), ...
    double(Tf.Ru1), double(Tf.Ru2), double(Tf.Ru3), ...
    double(Tf.Rdu1), double(Tf.Rdu2), double(Tf.Rdu3), ...
    noisy_Jtrack, noisy_JTV, ...
    noisy_settle_x1_h_c1, noisy_settle_x2_h_c1, noisy_settle_x3_h_c1, ...
    noisy_settle_x1_h_c2, noisy_settle_x2_h_c2, noisy_settle_x3_h_c2, ...
    noisy_settle_case_h_c1, noisy_settle_case_h_c2, ...
    noisy_IAE_x1_c1, noisy_IAE_x2_c1, noisy_IAE_x3_c1, ...
    noisy_IAE_x1_c2, noisy_IAE_x2_c2, noisy_IAE_x3_c2, ...
    noisy_IAE_case_c1, noisy_IAE_case_c2, ...
    noiseless_Jtrack, noiseless_JTV, ...
    noiseless_settle_x1_h_c1, noiseless_settle_x2_h_c1, noiseless_settle_x3_h_c1, ...
    noiseless_settle_x1_h_c2, noiseless_settle_x2_h_c2, noiseless_settle_x3_h_c2, ...
    noiseless_settle_case_h_c1, noiseless_settle_case_h_c2, ...
    noiseless_IAE_x1_c1, noiseless_IAE_x2_c1, noiseless_IAE_x3_c1, ...
    noiseless_IAE_x1_c2, noiseless_IAE_x2_c2, noiseless_IAE_x3_c2, ...
    noiseless_IAE_case_c1, noiseless_IAE_case_c2, ...
    'VariableNames', { ...
    'run_key','timestamp','Np','Nc', ...
    'Q1','Q2','Q3','Ru1','Ru2','Ru3','Rdu1','Rdu2','Rdu3', ...
    'noisy_Jtrack','noisy_JTV', ...
    'noisy_settle_x1_h_c1','noisy_settle_x2_h_c1','noisy_settle_x3_h_c1', ...
    'noisy_settle_x1_h_c2','noisy_settle_x2_h_c2','noisy_settle_x3_h_c2', ...
    'noisy_settle_case_h_c1','noisy_settle_case_h_c2', ...
    'noisy_IAE_x1_c1','noisy_IAE_x2_c1','noisy_IAE_x3_c1', ...
    'noisy_IAE_x1_c2','noisy_IAE_x2_c2','noisy_IAE_x3_c2', ...
    'noisy_IAE_case_c1','noisy_IAE_case_c2', ...
    'noiseless_Jtrack','noiseless_JTV', ...
    'noiseless_settle_x1_h_c1','noiseless_settle_x2_h_c1','noiseless_settle_x3_h_c1', ...
    'noiseless_settle_x1_h_c2','noiseless_settle_x2_h_c2','noiseless_settle_x3_h_c2', ...
    'noiseless_settle_case_h_c1','noiseless_settle_case_h_c2', ...
    'noiseless_IAE_x1_c1','noiseless_IAE_x2_c1','noiseless_IAE_x3_c1', ...
    'noiseless_IAE_x1_c2','noiseless_IAE_x2_c2','noiseless_IAE_x3_c2', ...
    'noiseless_IAE_case_c1','noiseless_IAE_case_c2'});

Treport = sortrows(Treport, 'noisy_JTV', 'descend');
fprintf('\nFinal Pareto frontier f=1 controller table (sorted by decreasing noisy J_TV):\n');
disp(Treport);

% Controllers simultaneously in top-3 (lowest) for noisy J_track and all
% individual noisy settling/IAE metrics.
noisyJtrackSel = double(Treport.noisy_Jtrack);
noisyJtrackSel(~isfinite(noisyJtrackSel)) = inf;

settleColNames = { ...
    'noisy_settle_x1_h_c1','noisy_settle_x2_h_c1','noisy_settle_x3_h_c1', ...
    'noisy_settle_x1_h_c2','noisy_settle_x2_h_c2','noisy_settle_x3_h_c2'};
iaeColNames = { ...
    'noisy_IAE_x1_c1','noisy_IAE_x2_c1','noisy_IAE_x3_c1', ...
    'noisy_IAE_x1_c2','noisy_IAE_x2_c2','noisy_IAE_x3_c2'};

nTop = min(3, height(Treport));
[~, ordJtrackSel] = sort(noisyJtrackSel, 'ascend');
topJtrackMask = false(height(Treport), 1);
topJtrackMask(ordJtrackSel(1:nTop)) = true;
topSettleAllMask = true(height(Treport), 1);
topIAEAllMask = true(height(Treport), 1);

for c = 1:numel(settleColNames)
    v = double(Treport.(settleColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, 'ascend');
    m = false(height(Treport), 1);
    m(ord(1:nTop)) = true;
    topSettleAllMask = topSettleAllMask & m;
end
for c = 1:numel(iaeColNames)
    v = double(Treport.(iaeColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, 'ascend');
    m = false(height(Treport), 1);
    m(ord(1:nTop)) = true;
    topIAEAllMask = topIAEAllMask & m;
end
topAllMask = topJtrackMask & topSettleAllMask & topIAEAllMask;

showCols = [{'run_key','timestamp','noisy_Jtrack'}, settleColNames, iaeColNames];
TtopAll = Treport(topAllMask, showCols);
if ~isempty(TtopAll)
    TtopAll = sortrows(TtopAll, 'noisy_Jtrack', 'ascend');
end
fprintf('Controllers simultaneously top-3 lowest in noisy J_track and all individual noisy settling/IAE metrics (NaN treated as Inf):\n');
disp(TtopAll);

[benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);
Tf1Bench = table();
if benchFound
    [nSuperiorF1, isSuperiorF1] = count_benchmark_strict_superior( ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, benchJtrack, benchJTV);
    fprintf('Final Pareto f=1 controllers strictly better than noisy benchmark (J_track and J_TV both lower): %d/%d\n', ...
        nSuperiorF1, height(Treport));
    print_benchmark_ratio_stats('Final Pareto f=1', ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, isSuperiorF1, benchJtrack, benchJTV);
    if nSuperiorF1 > 0
        fprintf('Timestamps strictly better than benchmark:\n');
        disp(Treport.timestamp(isSuperiorF1));
    end
    Tf1Bench = build_benchmark_comparison_table( ...
        double(Treport.noisy_Jtrack), double(Treport.noisy_JTV), benchJtrack, benchJTV);
    fprintf('Final Pareto f=1 benchmark comparison table:\n');
    disp(Tf1Bench);
else
    warning('Benchmark superiority verification skipped: noisy benchmark not available.');
end

if ~isfolder(numericalFolder)
    mkdir(numericalFolder);
end
writetable(Treport, fullfile(numericalFolder, 'final_pareto_frontier_f1_noisy_noiseless_metrics.csv'));

outTxtDir = fullfile(projectRoot, 'results', 'txt results');
if ~isfolder(outTxtDir)
    mkdir(outTxtDir);
end
if benchFound
    benchTxtPath = fullfile(outTxtDir, 'final_pareto_f1_benchmark_comparison_table.txt');
    writetable(Tf1Bench, benchTxtPath, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('Saved: %s\n', benchTxtPath);
end
end


function M = compute_out_metrics_by_path(matPath, settlingTol)
%COMPUTE_OUT_METRICS_BY_PATH Load one out_full MAT and compute objective/settling/IAE metrics.
M = struct();
M.Jtrack = nan;
M.JTV = nan;
M.settle_h = nan(2, 3);
M.settle_case_h = nan(2, 1);
M.IAE = nan(2, 3);
M.IAE_case = nan(2, 1);
if ~isfile(matPath)
    return
end
S = load(matPath, 'out');
if ~isfield(S, 'out')
    return
end
out = S.out;
if isfield(out, 'SSE'), M.Jtrack = double(out.SSE); end
if isfield(out, 'SSdU'), M.JTV = double(out.SSdU); end
if ~isfield(out, 'case') || isempty(out.case)
    return
end
nCase = min(numel(out.case), 2);
for c = 1:nCase
    caseData = out.case(c);
    [settle_h, iae] = summarize_case_metrics_simple(caseData, settlingTol);
    nState = min(numel(settle_h), 3);
    M.settle_h(c, 1:nState) = settle_h(1:nState);
    M.IAE(c, 1:nState) = iae(1:nState);
    if any(isfinite(settle_h))
        M.settle_case_h(c) = max(settle_h(isfinite(settle_h)));
    else
        M.settle_case_h(c) = nan;
    end
    M.IAE_case(c) = sum(iae, 'omitnan');
end
end


function [settlingTimes_h, IAEByState] = summarize_case_metrics_simple(caseStruct, settlingTol)
%SUMMARIZE_CASE_METRICS_SIMPLE Compute settling times and IAE from one case struct.
if ~isfield(caseStruct, 'Y') || ~isfield(caseStruct, 'Ysp')
    settlingTimes_h = nan(1, 3);
    IAEByState = nan(1, 3);
    return
end
Y = double(caseStruct.Y);
Ysp = double(caseStruct.Ysp);
nState = min(size(Y, 2), size(Ysp, 2));
Y = Y(:, 1:nState);
Ysp = Ysp(:, 1:nState);
if isfield(caseStruct, 'dt')
    dt = double(caseStruct.dt);
elseif isfield(caseStruct, 'tf')
    dt = double(caseStruct.tf) / max(size(Y, 1) - 1, 1);
else
    dt = 1/60;
end
t = (0:size(Y, 1) - 1).' * dt;
settlingTimes_h = nan(1, nState);
IAEByState = sum(abs(Y - Ysp), 1) * dt;
epsRef = 1e-9;
for i = 1:nState
    refVal = Ysp(end, i);
    relErr = abs(Y(:, i) - Ysp(:, i)) / max(abs(refVal), epsRef);
    if relErr(end) > settlingTol
        settlingTimes_h(i) = nan;
        continue
    end
    settleIdx = find_settling_index(relErr, settlingTol);
    if isempty(settleIdx)
        settlingTimes_h(i) = nan;
    else
        settlingTimes_h(i) = t(settleIdx);
    end
end
if nState < 3
    settlingTimes_h(1, end+1:3) = nan;
    IAEByState(1, end+1:3) = nan;
end
end


function idx = find_settling_index(relErr, tol)
%FIND_SETTLING_INDEX First index where relative error stays below tol thereafter.
n = numel(relErr);
idx = [];
for k = 1:n
    if all(relErr(k:end) <= tol)
        idx = k;
        return
    end
end
end
