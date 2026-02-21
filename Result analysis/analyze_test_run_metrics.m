%% ANALYZE_TEST_RUN_METRICS
% Offline diagnostics for NMPC replay files in results/test_run.
% No live simulation is executed.
%
% Settling-time definition (per state i):
%   e_rel_i(k) = |y_i(k) - r_i| / max(|r_i|, eps_ref),  eps_ref = 1e-9
%   t_s,i is the earliest time where e_rel_i(t) <= tol for ALL future t.
%   If final relative error is > tol, settling time is NaN.
%
% Figure-generation method summary
% - Load full-horizon replay MAT files from results/test_run.
% - Compute settling/error/IAE/input-variation metrics per case.
% - Build controller-level summaries (case 1 + case 2 joins).
% - Plot settling/final-error boxplots and final-Pareto `N_p` distributions.
% - Export CSV and text summaries used by downstream Python figure scripts.
%
clear; close all; clc;
scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
addpath(genpath(fullfile(projectRoot, "dependencies")));
NATURE_COLOR = nature_methods_colors();

% User-level configuration for input folders and analysis scope.
cfg = struct();
cfg.project_root = projectRoot;
cfg.results_root = fullfile(cfg.project_root, "results");
cfg.test_run_root = fullfile(cfg.results_root, "test_run");
cfg.analysis_dir = scriptDir;
cfg.run_folders = [
    struct("label","Case 1","subfolder","run1_full_f1_no_noise");
    struct("label","Case 2","subfolder","run2_full_f1_no_noise")
];
cfg.final_pareto_timestamps = [
    "20260131_151035";
    "20260201_111557";
    "20260201_111928";
    "20260201_192807";
    "20260201_223337";
    "20260201_232106";
    "20260210_151703";
    "20260210_171107";
    "20260210_180826";
    "20260211_122653";
    "20260211_134235"
];
cfg.settling_tol = 0.05;
cfg.optim_mean_sim_time_h = 8.09603;

diagnostics = run_full_run_diagnostics(cfg);
% Text summary first, then figures.
print_full_run_report(diagnostics, cfg);
write_numerical_results(diagnostics, cfg.settling_tol, cfg.optim_mean_sim_time_h, cfg.results_root);
export_boxplot_data_csv(diagnostics, cfg.analysis_dir);
plot_summary_boxplots(diagnostics, cfg.optim_mean_sim_time_h, NATURE_COLOR.Orange);
plot_p_distribution_by_run_pareto(diagnostics, cfg);

%% FUNCTIONS
function diagnostics = run_full_run_diagnostics(cfg)
%RUN_FULL_RUN_DIAGNOSTICS Aggregate case/controller metrics from replay MAT files.
arguments
    cfg struct
end

caseRecordList = [];
controllerRecordList = [];
matFileCount = 0;

% Traverse all configured replay folders (run1/run2).
for k = 1:numel(cfg.run_folders)
    runName = cfg.run_folders(k).label;
    runFolder = fullfile(cfg.test_run_root, cfg.run_folders(k).subfolder);
    if ~isfolder(runFolder)
        warning("Run folder missing: %s", runFolder);
        continue
    end

    % One replay artifact per timestamp/controller.
    matFiles = dir(fullfile(runFolder, "out_full_*.mat"));
    matFileCount = matFileCount + numel(matFiles);

    for m = 1:numel(matFiles)
        matFilePath = fullfile(runFolder, matFiles(m).name);
        loaded = load(matFilePath, "ts", "out");
        if ~isfield(loaded, "ts") || ~isfield(loaded, "out") || ~isfield(loaded.out, "case")
            warning("Skipping %s (missing ts/out/case).", matFilePath);
            continue
        end

        timestamp = string(loaded.ts);
        outputData = loaded.out;

        % Controller-level metadata (shared by both simulated cases).
        controllerRecord = struct();
        controllerRecord.run_label = string(runName);
        controllerRecord.timestamp = timestamp;
        controllerRecord.SSE = get_struct_field_or_nan(outputData, "SSE");
        controllerRecord.SSdU = get_struct_field_or_nan(outputData, "SSdU");
        controllerRecord.J = controllerRecord.SSE + 1e4 * controllerRecord.SSdU;
        controllerRecord.p = NaN;
        controllerRecord.m = NaN;
        controllerRecord.Q_diag = {NaN};
        controllerRecord.R1_diag = {NaN};
        controllerRecord.R2_diag = {NaN};
        [decodedTheta, thetaOk] = decode_theta_weights(outputData);
        if thetaOk
            controllerRecord.p = decodedTheta.p;
            controllerRecord.m = decodedTheta.m;
            controllerRecord.Q_diag = {decodedTheta.Q_diag};
            controllerRecord.R1_diag = {decodedTheta.R1_diag};
            controllerRecord.R2_diag = {decodedTheta.R2_diag};
        end
        controllerRecordList = [controllerRecordList; controllerRecord]; %#ok<AGROW>

        % Case-level metrics (two initial-condition scenarios per controller).
        for c = 1:numel(outputData.case)
            caseData = outputData.case(c);
            caseMetrics = summarize_case_metrics(outputData, caseData, cfg.settling_tol);

            caseRecord = struct();
            caseRecord.run_label = string(runName);
            caseRecord.timestamp = timestamp;
            caseRecord.case_id = caseData.case_id;
            caseRecord.tf_h = get_struct_field_or_nan(outputData, "tf");
            caseRecord.N_steps = size(caseData.Y, 1);
            caseRecord.SSE = controllerRecord.SSE;
            caseRecord.SSdU = controllerRecord.SSdU;
            caseRecord.J = controllerRecord.J;
            caseRecord.p = controllerRecord.p;
            caseRecord.m = controllerRecord.m;

            for s = 1:numel(caseMetrics.settling_times_h)
                caseRecord.(sprintf("settle_x%d_h", s)) = caseMetrics.settling_times_h(s);
                caseRecord.(sprintf("final_err_x%d_pct", s)) = caseMetrics.final_error_pct(s);
                caseRecord.(sprintf("peak_err_x%d_pct", s)) = caseMetrics.peak_error_pct(s);
                caseRecord.(sprintf("IAE_x%d", s)) = caseMetrics.IAE_by_state(s);
            end
            for u = 1:numel(caseMetrics.total_input_variation)
                caseRecord.(sprintf("tot_dU_u%d", u)) = caseMetrics.total_input_variation(u);
            end

            caseRecordList = [caseRecordList; caseRecord]; %#ok<AGROW>
        end
    end
end

% Convert accumulated structs to tables for downstream grouping/filtering.
if isempty(caseRecordList)
    perCaseTable = table();
else
    perCaseTable = struct2table(caseRecordList);
end
if isempty(controllerRecordList)
    perControllerTable = table();
else
    perControllerTable = struct2table(controllerRecordList);
    [~, uniqueIdx] = unique(perControllerTable(:, ["run_label","timestamp"]), "rows", "stable");
    perControllerTable = perControllerTable(uniqueIdx, :);
end

diagnostics = struct();
diagnostics.per_case = perCaseTable;
diagnostics.per_controller = perControllerTable;
diagnostics.mat_file_count = matFileCount;

% Early return path when no usable replay files are found.
if isempty(perCaseTable)
    diagnostics.pareto_cases = table();
    diagnostics.pareto_controllers = table();
    diagnostics.missing_pareto_timestamps = cfg.final_pareto_timestamps;
    return
end

diagnostics.pareto_cases = perCaseTable(ismember(perCaseTable.timestamp, cfg.final_pareto_timestamps), :);
diagnostics.pareto_controllers = perControllerTable(ismember(perControllerTable.timestamp, cfg.final_pareto_timestamps), :);
foundTimestamps = unique(perCaseTable.timestamp);
diagnostics.missing_pareto_timestamps = setdiff(cfg.final_pareto_timestamps, foundTimestamps);
end


function x = get_struct_field_or_nan(S, fieldName)
%GET_STRUCT_FIELD_OR_NAN Read optional struct field and return NaN if absent.
% Defensive accessor for optional fields in saved replay structs.
if isfield(S, fieldName)
    x = double(S.(fieldName));
else
    x = NaN;
end
end


function [decodedTheta, isValid] = decode_theta_weights(outputStruct)
%DECODE_THETA_WEIGHTS Decode theta vector into N_p/N_c and weight diagonals.
% Decode theta = [f, theta_p, theta_m, q(1:nx), r1(1:nu), r2(1:nu)].
decodedTheta = struct("p", NaN, "m", NaN, "Q_diag", [], "R1_diag", [], "R2_diag", []);
isValid = false;
if ~isfield(outputStruct, "theta")
    return
end

theta = double(outputStruct.theta(:).');
if numel(theta) < 4
    return
end

    nx = size(outputStruct.cfg.Q, 1);
    nu = size(outputStruct.cfg.Ru, 1);


expectedLength = 1 + 2 + nx + nu + nu;
if expectedLength ~= numel(theta)
    return
end

f_unused = theta(1); %#ok<NASGU>
thetaP = theta(2);
thetaM = theta(3);
qExp = theta(4:3+nx);
r1Exp = theta(4+nx:3+nx+nu);
r2Exp = theta(4+nx+nu:3+nx+2*nu);

decodedTheta.p = thetaP + (thetaM + 1);
decodedTheta.m = thetaM + 1;
decodedTheta.Q_diag = 10.^qExp;
decodedTheta.R1_diag = 10.^r1Exp;
decodedTheta.R2_diag = 10.^r2Exp;
isValid = true;
end


function metrics = summarize_case_metrics(outStruct, caseStruct, settlingTol)
%SUMMARIZE_CASE_METRICS Compute settling/error/input-variation metrics per case.
arguments
    outStruct struct
    caseStruct struct
    settlingTol (1,1) double {mustBePositive} = 0.02
end

stateTrajectory = double(caseStruct.Y);
setpointTrajectory = double(caseStruct.Ysp);
inputTrajectory = double(caseStruct.U);

if isfield(caseStruct, "dt")
    sampleTimeHours = double(caseStruct.dt);
elseif isfield(outStruct, "dt")
    sampleTimeHours = double(outStruct.dt);
else
    error("Missing dt in case/out struct.");
end
timeHours = (0:size(stateTrajectory,1)-1).' * sampleTimeHours;

[settlingTimes_h, finalErrPct, peakErrPct] = compute_settling_times( ...
    stateTrajectory, setpointTrajectory, timeHours, settlingTol);
totalInputVariation = compute_total_input_variation(inputTrajectory);
[IAEByState, ~] = compute_integral_abs_error(stateTrajectory, setpointTrajectory, sampleTimeHours);

metrics = struct();
metrics.settling_times_h = settlingTimes_h;
metrics.final_error_pct = finalErrPct;
metrics.peak_error_pct = peakErrPct;
metrics.IAE_by_state = IAEByState;
metrics.total_input_variation = totalInputVariation;
end


function [settlingTimes_h, finalErrPct, peakErrPct] = compute_settling_times(stateTrajectory, setpointTrajectory, timeHours, tol)
%COMPUTE_SETTLING_TIMES Compute per-state settling, final error, and peak error.
arguments
    stateTrajectory double
    setpointTrajectory double
    timeHours double
    tol (1,1) double {mustBeGreaterThan(tol,0)} = 0.02
end

numStates = size(stateTrajectory, 2);
settlingTimes_h = NaN(1, numStates);
finalErrPct = NaN(1, numStates);
peakErrPct = NaN(1, numStates);
eps_ref = 1e-9;

for stateIdx = 1:numStates
    refValue = setpointTrajectory(end, stateIdx);
    normalizer = max(abs(refValue), eps_ref);
    relError = abs(stateTrajectory(:, stateIdx) - refValue) ./ normalizer;

    finalErrPct(stateIdx) = 100 * relError(end);
    peakErrPct(stateIdx) = 100 * max(relError);

    if finalErrPct(stateIdx) > 100 * tol
        settlingTimes_h(stateIdx) = NaN;
        continue
    end

    % First time index from which all future samples remain within tolerance.
    futureMaxRelError = flipud(cummax(flipud(relError)));
    firstSettledIdx = find(futureMaxRelError <= tol, 1, "first");
    if ~isempty(firstSettledIdx)
        settlingTimes_h(stateIdx) = timeHours(firstSettledIdx);
    end
end
end


function totalInputVariation = compute_total_input_variation(inputTrajectory)
%COMPUTE_TOTAL_INPUT_VARIATION Sum absolute input increments per actuator.
% Sum of absolute input increments per actuator.
if isempty(inputTrajectory)
    totalInputVariation = [];
    return
end
if size(inputTrajectory,1) < 2
    totalInputVariation = zeros(1, size(inputTrajectory,2));
    return
end
inputIncrements = diff(inputTrajectory, 1, 1);
totalInputVariation = sum(abs(inputIncrements), 1);
end


function [IAEByState, absError] = compute_integral_abs_error(stateTrajectory, setpointTrajectory, sampleTimeHours)
%COMPUTE_INTEGRAL_ABS_ERROR Compute IAE per state on uniform sample grid.
% IAE per state with rectangular integration on uniform sampling.
absError = abs(stateTrajectory - setpointTrajectory);
IAEByState = sum(absError, 1) * sampleTimeHours;
end


function print_full_run_report(diagnostics, cfg)
%PRINT_FULL_RUN_REPORT Emit human-readable diagnostics summary to console.
caseTable = diagnostics.per_case;
controllerTable = diagnostics.per_controller;
settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));

% Print full controller table (both cases) before summary, with no removals.
allControllersBothCases = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, false);
if isempty(allControllersBothCases)
    fprintf("Controllers table (both case 1 and case 2, no removals): none\n");
else
    fprintf("Controllers table (both case 1 and case 2, no removals):\n");
    allControllersBothCases = sort_table_by_sse(allControllersBothCases);
    allCols = ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag"];
    settleBothCols = allControllersBothCases.Properties.VariableNames(startsWith(allControllersBothCases.Properties.VariableNames, "settle_"));
    allCols = [allCols, settleBothCols];
    allCols = allCols(ismember(allCols, allControllersBothCases.Properties.VariableNames));
    disp(allControllersBothCases(:, allCols));
end

fprintf("=== NMPC Full-Horizon Replay Report ===\n");
fprintf("Run folders: %s\n", strjoin(string({cfg.run_folders.subfolder}), ", "));
fprintf("MAT files scanned: %d\n", diagnostics.mat_file_count);
fprintf("Cases: %d | Controllers: %d\n", height(caseTable), height(controllerTable));

if isempty(caseTable)
    fprintf("No case data found.\n");
    return
end

for i = 1:numel(settlingColumns)
    nanCount = nnz(isnan(caseTable.(settlingColumns{i})));
    fprintf("NaN settling times %s: %d/%d\n", settlingColumns{i}, nanCount, height(caseTable));
end

% Case-level NaN diagnostics for settling times.
if ~isempty(settlingColumns)
    fprintf("Median settling times (h):\n");
    settlingData = table2array(caseTable(:, settlingColumns));
    medianSettling = nanmedian(settlingData, 1);
    for i = 1:numel(settlingColumns)
        fprintf("  %s: %.6g\n", settlingColumns{i}, medianSettling(i));
    end

    timeLimitH = cfg.optim_mean_sim_time_h;
    fprintf("Fraction settled before %.6g h:\n", timeLimitH);
    for i = 1:numel(settlingColumns)
        settledBefore = isfinite(caseTable.(settlingColumns{i})) & (caseTable.(settlingColumns{i}) <= timeLimitH);
        frac = nnz(settledBefore) / height(caseTable);
        fprintf("  %s: %.2f%% (%d/%d)\n", settlingColumns{i}, 100*frac, nnz(settledBefore), height(caseTable));
    end
    settledAllBefore = all(isfinite(settlingData) & (settlingData <= timeLimitH), 2);
    fracAll = nnz(settledAllBefore) / height(caseTable);
    fprintf("  all states: %.2f%% (%d/%d)\n", 100*fracAll, nnz(settledAllBefore), height(caseTable));
end

if ~isempty(diagnostics.missing_pareto_timestamps)
    fprintf("Missing configured Pareto timestamps: %s\n", strjoin(diagnostics.missing_pareto_timestamps, ", "));
end

if isempty(controllerTable)
    return
end

% Controller-level IAE table: one row per controller, with case-1/case-2 columns.
iaeColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "IAE_x"));
if ~isempty(iaeColumns)
    iaeCase1 = caseTable(caseTable.case_id == 1, ["run_label","timestamp",iaeColumns]);
    iaeCase2 = caseTable(caseTable.case_id == 2, ["run_label","timestamp",iaeColumns]);
    if ~isempty(iaeCase1) && ~isempty(iaeCase2)
        iaeCase1 = renamevars(iaeCase1, iaeColumns, iaeColumns + "_c1");
        iaeCase2 = renamevars(iaeCase2, iaeColumns, iaeColumns + "_c2");
        iaeControllerTable = innerjoin(iaeCase1, iaeCase2, "Keys", ["run_label","timestamp"]);
        iaeControllerTable = innerjoin( ...
            iaeControllerTable, ...
            controllerTable(:, ["run_label","timestamp","SSE"]), ...
            "Keys", ["run_label","timestamp"]);
        iaeControllerTable = sort_table_by_sse(iaeControllerTable);

        leftCols = ["timestamp","SSE","run_label"];
        iaeColsBoth = iaeControllerTable.Properties.VariableNames(startsWith(iaeControllerTable.Properties.VariableNames, "IAE_"));
        reportCols = [leftCols, iaeColsBoth];
        reportCols = reportCols(ismember(reportCols, iaeControllerTable.Properties.VariableNames));

        fprintf("Integral absolute error table (one row per controller):\n");
        disp(iaeControllerTable(:, reportCols));
    else
        fprintf("Integral absolute error table (one row per controller): none\n");
    end
end

% Controller-level listing: require BOTH cases and no NaN settling times in either case.
if ~isempty(settlingColumns)
    bothCasesTable = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, true);
    if isempty(bothCasesTable)
        fprintf("Controllers table (both cases with complete settling times): none\n");
    else
        bothCasesTable = sort_table_by_sse(bothCasesTable);
        settlingBothCols = bothCasesTable.Properties.VariableNames(startsWith(bothCasesTable.Properties.VariableNames, "settle_"));
        fprintf("Controllers table (both case 1 and case 2, excluding any NaN settling time):\n");
        reportCols = ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag",settlingBothCols];
        reportCols = reportCols(ismember(reportCols, bothCasesTable.Properties.VariableNames));
        disp(bothCasesTable(:, reportCols));
    end
end
end


function T = sort_table_by_sse(T)
%SORT_TABLE_BY_SSE Sort table rows by SSE ascending when column exists.
% Standard print ordering for all tables in this script.
if isempty(T) || ~ismember("SSE", string(T.Properties.VariableNames))
    return
end
T = sortrows(T, "SSE", "ascend");
end


function controllerCaseTable = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, removeNaNSettling)
%BUILD_CONTROLLER_BOTH_CASE_TABLE Join case-1 and case-2 metrics per controller.
% Build one row per controller by joining case 1 and case 2 side-by-side.
if isempty(caseTable) || isempty(controllerTable) || isempty(settlingColumns)
    controllerCaseTable = table();
    return
end

case1Table = caseTable(caseTable.case_id == 1, ["run_label","timestamp",settlingColumns]);
case2Table = caseTable(caseTable.case_id == 2, ["run_label","timestamp",settlingColumns]);
if isempty(case1Table) || isempty(case2Table)
    controllerCaseTable = table();
    return
end

% Suffix columns by case id to make the merged table explicit.
case1Table = renamevars(case1Table, settlingColumns, settlingColumns + "_c1");
case2Table = renamevars(case2Table, settlingColumns, settlingColumns + "_c2");

controllerCaseTable = innerjoin(case1Table, case2Table, "Keys", ["run_label","timestamp"]);
controllerCaseTable = innerjoin( ...
    controllerCaseTable, ...
    controllerTable(:, ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag"]), ...
    "Keys", ["run_label","timestamp"]);

if removeNaNSettling
    % Optional strict filter: drop controllers with any NaN settling metric.
    settlingBothCols = controllerCaseTable.Properties.VariableNames(startsWith(controllerCaseTable.Properties.VariableNames, "settle_"));
    validRowsMask = ~any(isnan(table2array(controllerCaseTable(:, settlingBothCols))), 2);
    controllerCaseTable = controllerCaseTable(validRowsMask, :);
end

if ~isempty(controllerCaseTable)
    controllerCaseTable = sortrows(controllerCaseTable, ["SSE","SSdU"], ["ascend","ascend"]);
end
end


function plot_summary_boxplots(diagnostics, settlingTimeRefH, guidelineColor)
%PLOT_SUMMARY_BOXPLOTS Plot settling-time and final-error distribution summaries.
caseTable = diagnostics.per_case;
if isempty(caseTable)
    return
end

% Compact diagnostic summary focused on settling time and final error.
figure("Color","w","Name","test_run_metric_summary");

settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));
if ~isempty(settlingColumns)
    subplot(1, 2, 1);
    stateLabels = get_state_display_labels(numel(settlingColumns));
    boxplot(table2array(caseTable(:, settlingColumns)), "Labels", stateLabels);
    yline(settlingTimeRefH, "-", "LineWidth", 1.2, "Color", guidelineColor);
    xlabel("State");
    ylabel("Settling time (h)");
    title("Settling Time");
    grid off; box off;
end

subplot(1, 2, 2);
errorColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "final_err_x"));
if ~isempty(errorColumns)
    errorLabels = get_state_display_labels(numel(errorColumns));
    boxplot(table2array(caseTable(:, errorColumns)), "Labels", errorLabels);
else
    boxplot(zeros(height(caseTable),1), "Labels", {'n/a'});
end
xlabel("State");
ylabel("Final relative error (%)");
title("Final Relative Error (%)");
grid off; box off;
end


function plot_p_distribution_by_run_pareto(diagnostics, cfg)
%PLOT_P_DISTRIBUTION_BY_RUN_PARETO Plot N_p distribution for final Pareto set.
caseTable = diagnostics.per_case;
controllerTable = diagnostics.per_controller;
if isempty(caseTable) || isempty(controllerTable)
    return
end

% Same run comparison as above, but restricted to configured final Pareto timestamps.
settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));
allControllersBothCases = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, false);
if isempty(allControllersBothCases) || ~ismember("p", string(allControllersBothCases.Properties.VariableNames))
    fprintf("Run comparison plot (p, Pareto): insufficient data.\n");
    return
end

isPareto = ismember(allControllersBothCases.timestamp, cfg.final_pareto_timestamps);
plotTable = allControllersBothCases(isPareto, :);
hasValidP = isfinite(plotTable.p);
plotTable = plotTable(hasValidP, :);
if isempty(plotTable)
    fprintf("Run comparison plot (p, Pareto): no finite p values.\n");
    return
end

runGroup = categorical(plotTable.run_label);
pVals = plotTable.p;

figure("Color","w","Name","p_distribution_by_run_final_pareto");

boxplot(pVals, runGroup);
xlabel("Case");
ylabel("Prediction horizon $N_p$");
title("Final Pareto Controllers: $N_p$ by case");
grid off; box off;
runNames = categories(runGroup);

fprintf("Run-wise median p (final Pareto controllers):\n");
for i = 1:numel(runNames)
    mask = (runGroup == runNames{i});
    fprintf("  %s: %.4g\n", string(runNames{i}), median(pVals(mask), "omitnan"));
end
end


function labels = get_state_display_labels(nStates)
%GET_STATE_DISPLAY_LABELS Map state index to publication labels.
baseLabels = {'V (L)', 'X (g/L)', 'S (g/L)'};
if nStates <= numel(baseLabels)
    labels = baseLabels(1:nStates);
else
    extraLabels = arrayfun(@(k) sprintf("x%d", k), 4:nStates, "UniformOutput", false);
    labels = [baseLabels, extraLabels];
end
end


function export_boxplot_data_csv(diagnostics, analysisDir)
%EXPORT_BOXPLOT_DATA_CSV Write long-form CSVs consumed by Python plotting script.
caseTable = diagnostics.per_case;
if isempty(caseTable)
    return
end
if ~isfolder(analysisDir)
    mkdir(analysisDir);
end

settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));
errorColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "final_err_x"));

settlingLong = stack(caseTable(:, ["run_label","timestamp","case_id",settlingColumns]), ...
    settlingColumns, "NewDataVariableName", "value", "IndexVariableName", "source_col");
settlingLong.metric = repmat("settling_time_h", height(settlingLong), 1);
settlingLong.state = extractBetween(string(settlingLong.source_col), "settle_", "_h");
settlingLong.state_label = map_state_labels_from_id(settlingLong.state);
settlingLong = movevars(settlingLong, ["metric","state","state_label","value"], "After", "case_id");
settlingLong = removevars(settlingLong, "source_col");

errorLong = stack(caseTable(:, ["run_label","timestamp","case_id",errorColumns]), ...
    errorColumns, "NewDataVariableName", "value", "IndexVariableName", "source_col");
errorLong.metric = repmat("final_relative_error_pct", height(errorLong), 1);
errorLong.state = extractBetween(string(errorLong.source_col), "final_err_", "_pct");
errorLong.state_label = map_state_labels_from_id(errorLong.state);
errorLong = movevars(errorLong, ["metric","state","state_label","value"], "After", "case_id");
errorLong = removevars(errorLong, "source_col");

writetable(settlingLong, fullfile(analysisDir, "boxplot_settling_time_data.csv"));
writetable(errorLong, fullfile(analysisDir, "boxplot_final_error_data.csv"));
end


function stateLabels = map_state_labels_from_id(stateId)
%MAP_STATE_LABELS_FROM_ID Convert state IDs (x1/x2/x3/...) to display labels.
stateLabels = strings(numel(stateId), 1);
for i = 1:numel(stateId)
    switch string(stateId(i))
        case "x1"
            stateLabels(i) = "V (L)";
        case "x2"
            stateLabels(i) = "X (g/L)";
        case "x3"
            stateLabels(i) = "S (g/L)";
        otherwise
            stateLabels(i) = string(stateId(i));
    end
end
end


function write_numerical_results(diagnostics, settlingTol, settlingTimeRefH, resultsRoot)
%WRITE_NUMERICAL_RESULTS Persist computed diagnostics and summary statistics.
caseTable = diagnostics.per_case;
if isempty(caseTable)
    return
end

settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));
if isempty(settlingColumns)
    return
end

outDir = fullfile(resultsRoot, "numerical results");
if ~isfolder(outDir)
    mkdir(outDir);
end
outPath = fullfile(outDir, "analyze_test_run_metrics_summary.txt");
fid = fopen(outPath, "w");
if fid == -1
    warning("Unable to write numerical summary: %s", outPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for i = 1:numel(settlingColumns)
    col = settlingColumns{i};
    nanCount = nnz(isnan(caseTable.(col)));
    totalCount = height(caseTable);
    pct = 100 * nanCount / max(totalCount, 1);
    fprintf(fid, "NaN settling times %s: %d/%d\n", col, nanCount, totalCount);
    fprintf(fid, "Interpretation: %.2f%% of runs had final error > %.2f%% for %s.\n", ...
        pct, 100 * settlingTol, col);
end

fprintf(fid, "\nFraction settled before %.6g h:\n", settlingTimeRefH);
for i = 1:numel(settlingColumns)
    col = settlingColumns{i};
    settledBefore = isfinite(caseTable.(col)) & (caseTable.(col) <= settlingTimeRefH);
    frac = nnz(settledBefore) / height(caseTable);
    fprintf(fid, "  %s: %.2f%% (%d/%d)\n", col, 100*frac, nnz(settledBefore), height(caseTable));
end
settlingData = table2array(caseTable(:, settlingColumns));
settledAllBefore = all(isfinite(settlingData) & (settlingData <= settlingTimeRefH), 2);
fracAll = nnz(settledAllBefore) / height(caseTable);
fprintf(fid, "  all states: %.2f%% (%d/%d)\n", 100*fracAll, nnz(settledAllBefore), height(caseTable));
end
