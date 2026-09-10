function run_refined_frontier_change(E, caseNames, caseLabels, projectRoot, graphicsFolder, ...
    plotColors, caseMarkers, accentColor, fontSize)
%RUN_REFINED_FRONTIER_CHANGE Compare original (BO) and refined (z=1 re-eval) Pareto points.
% Guarded: skips with a warning when no refined result CSVs exist.
nCases = numel(E);
resultsRoot = fullfile(projectRoot, 'results');

refinedFiles = strings(nCases, 1);
for k = 1:nCases
    refinedFiles(k) = fullfile(resultsRoot, 'final_fidelity_same_noise', ...
        string(caseNames{k}) + "_full_f1_same_noise", 'results_full.csv');
end
present = arrayfun(@isfile, refinedFiles);
if ~any(present)
    warning('main_pareto:noRefined', ...
        ['Refined frontier section skipped: no z=1 re-evaluation results found.\n' ...
         'Expected e.g. %s'], refinedFiles(1));
    return
end
if ~all(present)
    warning('main_pareto:partialRefined', ...
        'Refined results missing for %d case(s); comparing the available ones only.', ...
        nnz(~present));
end

% Original pool: all BO evaluations, keyed case|timestamp.
T_orig = table();
for k = 1:nCases
    Tk = E{k};
    Tk.run_key   = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.run_label = repmat(caseLabels(k), height(Tk), 1);
    Tk.ts        = ts_string(Tk.timestamp);
    Tk.match_key = Tk.run_key + "|" + Tk.ts;
    Tk.z_eval    = double(Tk.z);
    T_orig = [T_orig; Tk(:, ["run_key", "run_label", "match_key", "ts", ...
        "iter", "SSE", "SSdU", "z_eval"])]; %#ok<AGROW>
end

% Refined pool from the re-evaluation CSVs.
T_ref = table();
for k = 1:nCases
    if ~present(k); continue; end
    Tk = readtable(refinedFiles(k), 'TextType', 'string');
    required = ["timestamp", "SSE", "SSdU"];
    for c = required
        if ~ismember(c, string(Tk.Properties.VariableNames))
            error('Missing column ''%s'' in %s', c, refinedFiles(k));
        end
    end
    Tk.run_key   = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.run_label = repmat(caseLabels(k), height(Tk), 1);
    Tk.ts        = string(Tk.timestamp);
    Tk.match_key = Tk.run_key + "|" + Tk.ts;
    T_ref = [T_ref; Tk(:, ["run_key", "run_label", "match_key", "ts", "SSE", "SSdU"])]; %#ok<AGROW>
end
if isempty(T_ref)
    warning('main_pareto:emptyRefined', 'Refined result files were found but contained no rows.');
    return
end

isParetoOrig = compute_pareto_mask(double(T_orig.SSE), double(T_orig.SSdU));
[T_ref, refFilterInfo] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig);
print_refined_filter_summary(refFilterInfo);
if isempty(T_ref)
    warning('main_pareto:noMatchedRefined', ...
        'No refined rows match original Pareto points; refined section skipped.');
    return
end
origParetoKeys = unique(string(T_orig.match_key(isParetoOrig)), 'stable');
missingRefined = setdiff(origParetoKeys, unique(string(T_ref.match_key), 'stable'), 'stable');
if ~isempty(missingRefined)
    warning('main_pareto:missingRefined', ...
        'Missing z=1 refined results for %d Pareto controller(s): %s', ...
        numel(missingRefined), strjoin(missingRefined, ', '));
end

T_jtv = compute_change_table(T_orig, T_ref, 'SSdU', 'JTV');
fprintf('\n=== J_TV change (original -> refined), sorted by |delta %%| descending ===\n');
disp(T_jtv);
T_jtrack = compute_change_table(T_orig, T_ref, 'SSE', 'Jtrack');
fprintf('\n=== J_track change (original -> refined), sorted by |delta %%| descending ===\n');
disp(T_jtrack);
print_top1_cfg(projectRoot, T_jtv);

% Promotion bookkeeping.
isParetoRef = compute_pareto_mask(double(T_ref.SSE), double(T_ref.SSdU));
[commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
[matchedRunKey, matchedTs] = split_match_key(commonKeys);

origParetoMatched = isParetoOrig(idxOrig);
refParetoMatched  = isParetoRef(idxRef);
zOrigMatched      = double(T_orig.z_eval(idxOrig));
promotedMask      = ~origParetoMatched & refParetoMatched;
promotedZlt1Mask  = promotedMask & isfinite(zOrigMatched) & (zOrigMatched < 1 - 1e-12);
promotedIdxOrig   = idxOrig(promotedZlt1Mask);
promotedIdxRef    = idxRef(promotedZlt1Mask);
promotedRunKey    = matchedRunKey(promotedZlt1Mask);
promotedTs        = matchedTs(promotedZlt1Mask);

% Per-case Pareto subsets of the original pool (for the base plot).
Ecase = cell(nCases, 1);
Tpcase = cell(nCases, 1);
for k = 1:nCases
    Tk = T_orig(string(T_orig.run_key) == string(caseNames{k}), :);
    Ecase{k} = Tk;
    Tpcase{k} = Tk(compute_pareto_mask(double(Tk.SSE), double(Tk.SSdU)), :);
end

NATURE_COLOR = nature_methods_colors();
colRefAll   = [0.60, 0.82, 0.98];          % re-evaluated reference cloud
colPromoted = plotColors(2, :);            % "promoted" points
benchmarkColor = [242, 133, 34] / 255;     % good_colors.m C.orange
[benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);

% Axis limits from the pooled original + refined data.
xLimR = padded_log_limits([double(T_orig.SSdU); double(T_ref.SSdU)], 1.25);
yLimR = padded_log_limits([double(T_orig.SSE);  double(T_ref.SSE)],  1.25);
if benchFound
    xLimR = expand_log_limits_to_include_point(xLimR, benchJTV, 1.10);
    yLimR = expand_log_limits_to_include_point(yLimR, benchJtrack, 1.10);
end

fig = figure('Color', 'w', 'Toolbar', 'none', 'Name', 'Pareto Frontier Change');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
set(fig, 'Position', [80, 80, 1400, 560]);

axL = nexttile; hold(axL, 'on');
[isParetoLeft, ~] = plot_combined_pareto_base(axL, Ecase, Tpcase, plotColors, caseMarkers, ...
    accentColor, xLimR, yLimR);
if benchFound
    scatter(axL, benchJTV, benchJtrack, 130, 's', ...
        'MarkerFaceColor', benchmarkColor, 'MarkerEdgeColor', 'none', 'LineWidth', 1.0);
end

axR = nexttile; hold(axR, 'on');
plot_combined_samples_no_guide(axR, Ecase, Tpcase, plotColors, caseMarkers);
scatter(axR, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
    'filled', 'MarkerFaceColor', colRefAll, 'MarkerEdgeColor', 'none');
TrefPareto = T_ref(isParetoRef, :);
scatter(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), 80, ...
    'd', 'MarkerFaceColor', NATURE_COLOR.ReddishPurple, ...
    'MarkerEdgeColor', NATURE_COLOR.ReddishPurple, 'LineWidth', 1.0);
if ~isempty(promotedIdxRef)
    scatter(axR, double(T_ref.SSdU(promotedIdxRef)), double(T_ref.SSE(promotedIdxRef)), 140, ...
        'p', 'MarkerFaceColor', colPromoted, 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
end
if benchFound
    scatter(axR, benchJTV, benchJtrack, 130, 's', ...
        'MarkerFaceColor', benchmarkColor, 'MarkerEdgeColor', 'none', 'LineWidth', 1.0);
end

if benchFound
    [nOrigSuperior, superiorMaskOrig] = count_benchmark_strict_superior( ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
    [nRefSuperior, superiorMaskRef] = count_benchmark_strict_superior( ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
    fprintf('Pareto controllers strictly better than noisy benchmark (J_track and J_TV both lower):\n');
    fprintf('  Original Pareto (BO phase): %d/%d\n', nOrigSuperior, nnz(isParetoOrig));
    print_benchmark_ratio_stats('Original Pareto (BO phase)', ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), ...
        superiorMaskOrig, benchJtrack, benchJTV);
    fprintf('  Refined Pareto: %d/%d\n', nRefSuperior, height(TrefPareto));
    print_benchmark_ratio_stats('Refined Pareto', ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), superiorMaskRef, benchJtrack, benchJTV);
    TorigBench = build_benchmark_comparison_table( ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
    fprintf('  Original Pareto benchmark comparison table:\n');
    disp(TorigBench);
    TrefBench = build_benchmark_comparison_table( ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
    fprintf('  Refined Pareto benchmark comparison table:\n');
    disp(TrefBench);
end

plot_pareto_continuum_line_only(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), ...
    NATURE_COLOR.ReddishPurple, xLimR, yLimR);
apply_combined_axes_style(axL, fontSize, xLimR, yLimR);
apply_combined_axes_style(axR, fontSize, xLimR, yLimR);

title(axL, '$\mathbf{a}$', 'Interpreter', 'latex');
title(axR, '$\mathbf{b}$', 'Interpreter', 'latex');
axL.TitleHorizontalAlignment = 'left';
axR.TitleHorizontalAlignment = 'left';

outStem = fullfile(graphicsFolder, 'refined_frontier_change');
exportgraphics(fig, outStem + ".png", 'Resolution', 300);
exportgraphics(fig, outStem + ".pdf", 'ContentType', 'vector');

outTxtDir = fullfile(projectRoot, 'results', 'txt results');
if ~isfolder(outTxtDir); mkdir(outTxtDir); end
reportPath = fullfile(outTxtDir, 'refined_promoted_frontier_z_lt_1.txt');
write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, ...
    promotedIdxOrig, promotedIdxRef);

fprintf('Saved: %s\n', outStem + ".png");
fprintf('Saved: %s\n', outStem + ".pdf");
fprintf('Saved: %s\n', reportPath);
fprintf('J_TV rows compared: %d\n', height(T_jtv));
fprintf('Matched points (case + timestamp): %d\n', numel(idxOrig));
fprintf('Left-panel frontier pool size (all BO points): %d\n', height(T_orig));
fprintf('Left-panel frontier points (combined Pareto): %d\n', nnz(isParetoLeft));
fprintf('Original Pareto points (all BO rows): %d\n', nnz(isParetoOrig));
fprintf('Refined Pareto points (all refined rows): %d\n', nnz(isParetoRef));
fprintf('Promoted to Pareto with z < 1: %d\n', numel(promotedIdxRef));
end


function [TrefKeep, info] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig)
%KEEP_REFINED_FROM_ORIGINAL_PARETO Keep only refined points linked to original Pareto.
refKeys        = string(T_ref.match_key);
origKeys       = string(T_orig.match_key);
origParetoKeys = string(T_orig.match_key(isParetoOrig));

inOrig       = ismember(refKeys, origKeys);
inOrigPareto = ismember(refKeys, origParetoKeys);

info = struct();
info.n_ref_input            = height(T_ref);
info.n_drop_not_in_orig     = nnz(~inOrig);      % includes DOE-linked rows (not in evals)
info.n_drop_in_orig_not_pareto = nnz(inOrig & ~inOrigPareto);
TrefKeep = T_ref(inOrigPareto, :);
if ~isempty(TrefKeep)
    [~, uniqIdx] = unique(string(TrefKeep.match_key), 'stable');
    info.n_drop_duplicate_refined = height(TrefKeep) - numel(uniqIdx);
    TrefKeep = TrefKeep(uniqIdx, :);
else
    info.n_drop_duplicate_refined = 0;
end
info.n_keep_final = height(TrefKeep);
end


function print_refined_filter_summary(info)
%PRINT_REFINED_FILTER_SUMMARY Print keep/drop diagnostics.
fprintf('\n=== Refined sample eligibility check ===\n');
fprintf('Input refined rows: %d\n', info.n_ref_input);
fprintf('Dropped (not a BO evaluation: DOE-linked or unknown): %d\n', info.n_drop_not_in_orig);
fprintf('Dropped (BO evaluation but not Pareto): %d\n', info.n_drop_in_orig_not_pareto);
fprintf('Dropped duplicate refined rows: %d\n', info.n_drop_duplicate_refined);
fprintf('Kept refined rows (BO Pareto-linked): %d\n', info.n_keep_final);
end


function T = compute_change_table(T_orig, T_ref, col, label)
%COMPUTE_CHANGE_TABLE Build a delta table (original -> refined) for one cost column.
[commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
if isempty(commonKeys)
    error('No matching (case, timestamp) points between original and refined tables.');
end
[runKey, ts] = split_match_key(commonKeys);
T = table(runKey, ts, double(T_orig.(col)(idxOrig)), double(T_ref.(col)(idxRef)), ...
    'VariableNames', ["run_key", "timestamp", label + "_original", label + "_refined"]);
T.("delta_" + label) = T.(label + "_refined") - T.(label + "_original");
T.delta_pct = 100 * T.("delta_" + label) ./ max(abs(T.(label + "_original")), eps);
T.("abs_delta_" + label) = abs(T.("delta_" + label));
T.abs_delta_pct = abs(T.delta_pct);
T = sortrows(T, 'abs_delta_pct', 'descend');
end


function print_top1_cfg(projectRoot, T_jtv)
%PRINT_TOP1_CFG Print original/refined configs for the largest JTV change.
if isempty(T_jtv)
    return
end
runKey = string(T_jtv.run_key(1));
ts     = string(T_jtv.timestamp(1));
origMat = find_mat_for_timestamp(projectRoot, runKey, ts, true);
refMat  = find_mat_for_timestamp(projectRoot, runKey, ts, false);
if strlength(origMat) == 0 || strlength(refMat) == 0
    fprintf('\nTop-1 cfg print skipped (missing MAT files) for %s | %s.\n', runKey, ts);
    return
end
Sorig = load_cfg_snapshot(origMat);
Sref  = load_cfg_snapshot(refMat);
fprintf('\n=== Top-1 point by |delta J_TV|: %s | %s ===\n', runKey, ts);
fprintf('\n--- Original cfg (%s) ---\n', origMat);
print_cfg_struct(Sorig);
fprintf('\n--- Refined cfg (%s) ---\n', refMat);
print_cfg_struct(Sref);
end


function matPath = find_mat_for_timestamp(projectRoot, runKey, ts, isOriginal)
%FIND_MAT_FOR_TIMESTAMP Resolve original/refined MAT path for one controller.
if isOriginal
    cand = fullfile(projectRoot, 'results', runKey, "out_" + ts + ".mat");
else
    cand = fullfile(projectRoot, 'results', 'final_fidelity_same_noise', ...
        runKey + "_full_f1_same_noise", "out_full_" + ts + ".mat");
end
if isfile(cand)
    matPath = string(cand);
else
    matPath = "";
end
end


function S = load_cfg_snapshot(matPath)
%LOAD_CFG_SNAPSHOT Load only cfg-relevant MAT fields.
S = struct();
vars = string(who('-file', matPath));
if any(vars == "out")
    tmp = load(matPath, 'out');
    S.out = tmp.out;
end
if any(vars == "cfg_run")
    tmp = load(matPath, 'cfg_run');
    S.cfg_run = tmp.cfg_run;
end
end


function print_cfg_struct(S)
%PRINT_CFG_STRUCT Print compact configuration information.
if isfield(S, 'out') && isfield(S.out, 'cfg')
    cfg = S.out.cfg;
    if isfield(cfg, 'f'), fprintf('f = %.12g\n', cfg.f); end
    if isfield(cfg, 'm'), fprintf('m = %d\n', cfg.m); end
    if isfield(cfg, 'p'), fprintf('p = %d\n', cfg.p); end
    if isfield(cfg, 'Q'),   fprintf('Q =\n');   disp(cfg.Q);   end
    if isfield(cfg, 'Ru'),  fprintf('Ru =\n');  disp(cfg.Ru);  end
    if isfield(cfg, 'Rdu'), fprintf('Rdu =\n'); disp(cfg.Rdu); end
end
if isfield(S, 'out') && isfield(S.out, 'theta')
    fprintf('theta =\n');
    disp(S.out.theta);
end
if isfield(S, 'cfg_run')
    cfgRun = S.cfg_run;
    if isfield(cfgRun, 'sigma_y'), fprintf('sigma_y = '); disp(cfgRun.sigma_y); end
    if isfield(cfgRun, 'mode'), fprintf('cfg_run.mode = %s\n', string(cfgRun.mode)); end
    if isfield(cfgRun, 'source_root'), fprintf('cfg_run.source_root = %s\n', string(cfgRun.source_root)); end
    if isfield(cfgRun, 'output_root'), fprintf('cfg_run.output_root = %s\n', string(cfgRun.output_root)); end
    if isfield(cfgRun, 'results_csv'), fprintf('cfg_run.results_csv = %s\n', string(cfgRun.results_csv)); end
    if isfield(cfgRun, 'out_dir'), fprintf('cfg_run.out_dir = %s\n', string(cfgRun.out_dir)); end
    if isfield(cfgRun, 'NumWorkers'), fprintf('cfg_run.NumWorkers = %d\n', double(cfgRun.NumWorkers)); end
end
end


function [runKey, ts] = split_match_key(matchKey)
%SPLIT_MATCH_KEY Split "case|timestamp" identifiers.
n = numel(matchKey);
runKey = strings(n, 1);
ts = strings(n, 1);
for i = 1:n
    parts = split(string(matchKey(i)), '|');
    runKey(i) = parts(1);
    ts(i) = parts(2);
end
end


function write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, ...
    promotedIdxOrig, promotedIdxRef)
%WRITE_PROMOTED_REPORT Save promoted-point table for traceability.
fid = fopen(reportPath, 'w');
if fid < 0
    warning('Could not write promoted-point report: %s', reportPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Promoted to refined Pareto frontier with original z < 1\n');
fprintf(fid, 'Original pool is the BO phase only (DOE excluded structurally).\n');
fprintf(fid, 'Definition: original non-Pareto -> refined Pareto AND original z < 1\n');
fprintf(fid, 'Count: %d\n\n', numel(promotedTs));
fprintf(fid, 'case,timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n');
for i = 1:numel(promotedTs)
    io = promotedIdxOrig(i);
    ir = promotedIdxRef(i);
    fprintf(fid, '%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
        promotedRunKey(i), promotedTs(i), ...
        double(T_orig.SSE(io)), double(T_orig.SSdU(io)), double(T_orig.z_eval(io)), ...
        double(T_ref.SSE(ir)), double(T_ref.SSdU(ir)));
end
end


function limOut = expand_log_limits_to_include_point(limIn, v, padFactor)
%EXPAND_LOG_LIMITS_TO_INCLUDE_POINT Expand positive log-axis limits to include point v.
limOut = limIn;
if ~(isfinite(v) && v > 0)
    return
end
if nargin < 3 || ~isfinite(padFactor) || padFactor <= 1
    padFactor = 1.10;
end
if v < limOut(1)
    limOut(1) = v / padFactor;
end
if v > limOut(2)
    limOut(2) = v * padFactor;
end
end

%% ===================== GUARDED: BENCHMARK =====================
