%% ANALYZE_TEST_RUN_METRICS
% Offline diagnostics for NMPC replay files in results/test_run.
% No live simulation is executed.
%
% Settling-time definition (per state i):
%   e_rel_i(k) = |y_i(k) - r_i| / max(|r_i|, eps_ref),  eps_ref = 1e-9
%   t_s,i is the earliest time where e_rel_i(t) <= tol for ALL future t.
%   If final relative error is > tol, settling time is NaN.
%
% Observed stability criterion used in this script:
%   Case is "possibly stable" if the approximated value function is
%   nonincreasing along the trajectory: V(k+1) <= V(k) for all k.
%   Controller (timestamp) is stable if all its cases are stable.

clear; clc;

cfg = struct();
cfg.test_run_root = fullfile("results", "test_run");
cfg.run_folders = [
    struct("label","run1","subfolder","run1_full_f1_no_noise");
    struct("label","run2","subfolder","run2_full_f1_no_noise")
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

diagnostics = run_full_run_diagnostics(cfg);
print_full_run_report(diagnostics, cfg);
plot_summary_boxplots(diagnostics);

%% FUNCTIONS
function diagnostics = run_full_run_diagnostics(cfg)
arguments
    cfg struct
end

caseRecordList = [];
controllerRecordList = [];
matFileCount = 0;

for k = 1:numel(cfg.run_folders)
    runName = cfg.run_folders(k).label;
    runFolder = fullfile(cfg.test_run_root, cfg.run_folders(k).subfolder);
    if ~isfolder(runFolder)
        warning("Run folder missing: %s", runFolder);
        continue
    end

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
            end
            for u = 1:numel(caseMetrics.total_input_variation)
                caseRecord.(sprintf("tot_dU_u%d", u)) = caseMetrics.total_input_variation(u);
            end

            caseRecord.value_V1 = caseMetrics.value_V1;
            caseRecord.value_VN = caseMetrics.value_VN;
            caseRecord.stage_cost_sum = caseMetrics.stage_cost_sum;
            caseRecord.value_noninc_ratio = caseMetrics.value_noninc_ratio;
            caseRecord.value_case_stable = caseMetrics.value_case_stable;

            caseRecordList = [caseRecordList; caseRecord]; %#ok<AGROW>
        end
    end
end

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

if isempty(perCaseTable)
    diagnostics.pareto_cases = table();
    diagnostics.pareto_controllers = table();
    diagnostics.missing_pareto_timestamps = cfg.final_pareto_timestamps;
    diagnostics.controller_stability = table();
    diagnostics.stable_low_quartile = table();
    return
end

diagnostics.pareto_cases = perCaseTable(ismember(perCaseTable.timestamp, cfg.final_pareto_timestamps), :);
diagnostics.pareto_controllers = perControllerTable(ismember(perControllerTable.timestamp, cfg.final_pareto_timestamps), :);
foundTimestamps = unique(perCaseTable.timestamp);
diagnostics.missing_pareto_timestamps = setdiff(cfg.final_pareto_timestamps, foundTimestamps);

diagnostics.controller_stability = compute_controller_stability(perCaseTable);
diagnostics.stable_low_quartile = select_stable_low_quartile(diagnostics.controller_stability, perControllerTable);
end


function x = get_struct_field_or_nan(S, fieldName)
if isfield(S, fieldName)
    x = double(S.(fieldName));
else
    x = NaN;
end
end


function [decodedTheta, isValid] = decode_theta_weights(outputStruct)
decodedTheta = struct("p", NaN, "m", NaN, "Q_diag", [], "R1_diag", [], "R2_diag", []);
isValid = false;
if ~isfield(outputStruct, "theta")
    return
end

theta = double(outputStruct.theta(:).');
if numel(theta) < 4
    return
end

if isfield(outputStruct, "cfg") && isfield(outputStruct.cfg, "Q") && isfield(outputStruct.cfg, "Ru")
    nx = size(outputStruct.cfg.Q, 1);
    nu = size(outputStruct.cfg.Ru, 1);
else
    nx = floor((numel(theta) - 3) / 3);
    nu = nx;
end
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


function approxValue = compute_approx_value_function(stateTrajectory, setpointTrajectory, inputTrajectory, decodedTheta)
numSteps = size(stateTrajectory, 1);
trackingError = stateTrajectory - setpointTrajectory;

stateCost = sum((trackingError.^2) .* decodedTheta.Q_diag(:).', 2);
inputCost = sum((inputTrajectory.^2) .* decodedTheta.R1_diag(:).', 2);

inputIncrements = [zeros(1, size(inputTrajectory,2)); diff(inputTrajectory, 1, 1)];
deltaInputCost = sum((inputIncrements.^2) .* decodedTheta.R2_diag(:).', 2);

stageCost = stateCost + inputCost + deltaInputCost;

% Backward approximation with zero terminal value: V_{N+1} ~= 0.
% This implies V_k ~= sum_{j=k..N} stageCost(j). We do not validate this assumption here.
valueSequence = zeros(numSteps, 1);
runningTailCost = 0;
for k = numSteps:-1:1
    runningTailCost = runningTailCost + stageCost(k);
    valueSequence(k) = runningTailCost;
end

approxValue = struct();
approxValue.stage_cost = stageCost;
approxValue.V = valueSequence;
end


function metrics = summarize_case_metrics(outStruct, caseStruct, settlingTol)
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
[decodedTheta, thetaOk] = decode_theta_weights(outStruct);
if thetaOk
    approxValue = compute_approx_value_function(stateTrajectory, setpointTrajectory, inputTrajectory, decodedTheta);
    valueDiff = diff(approxValue.V);
    finiteValueDiff = valueDiff(isfinite(valueDiff));
    if isempty(finiteValueDiff)
        valueNonIncRatio = NaN;
        isValueCaseStable = false;
    else
        valueNonIncRatio = mean(finiteValueDiff <= 0);
        isValueCaseStable = all(finiteValueDiff <= 0);
    end
    valueV1 = approxValue.V(1);
    valueVN = approxValue.V(end);
    stageCostSum = sum(approxValue.stage_cost, "omitnan");
else
    valueV1 = NaN;
    valueVN = NaN;
    stageCostSum = NaN;
    valueNonIncRatio = NaN;
    isValueCaseStable = false;
end

metrics = struct();
metrics.settling_times_h = settlingTimes_h;
metrics.final_error_pct = finalErrPct;
metrics.peak_error_pct = peakErrPct;
metrics.total_input_variation = totalInputVariation;
metrics.value_V1 = valueV1;
metrics.value_VN = valueVN;
metrics.stage_cost_sum = stageCostSum;
metrics.value_noninc_ratio = valueNonIncRatio;
metrics.value_case_stable = isValueCaseStable;
end


function [settlingTimes_h, finalErrPct, peakErrPct] = compute_settling_times(stateTrajectory, setpointTrajectory, timeHours, tol)
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

    futureMaxRelError = flipud(cummax(flipud(relError)));
    firstSettledIdx = find(futureMaxRelError <= tol, 1, "first");
    if ~isempty(firstSettledIdx)
        settlingTimes_h(stateIdx) = timeHours(firstSettledIdx);
    end
end
end


function totalInputVariation = compute_total_input_variation(inputTrajectory)
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


function controllerStability = compute_controller_stability(caseTable)
grouped = groupsummary(caseTable, ["run_label","timestamp"], "all", "value_case_stable");
controllerStability = table();
controllerStability.run_label = grouped.run_label;
controllerStability.timestamp = grouped.timestamp;
controllerStability.case_count = grouped.GroupCount;
controllerStability.stable_case_count = grouped.sum_value_case_stable;
controllerStability.all_cases_stable = (grouped.sum_value_case_stable == grouped.GroupCount);
end


function stableLowQ = select_stable_low_quartile(controllerStability, controllerTable)
if isempty(controllerStability) || isempty(controllerTable)
    stableLowQ = table();
    return
end

validSSE = controllerTable.SSE(isfinite(controllerTable.SSE));
validSSdU = controllerTable.SSdU(isfinite(controllerTable.SSdU));
if isempty(validSSE) || isempty(validSSdU)
    stableLowQ = table();
    return
end
trackingCostQ1 = prctile(validSSE, 25);
inputUseCostQ1 = prctile(validSSdU, 25);

joined = innerjoin(controllerStability, controllerTable, "Keys", ["run_label","timestamp"]);
isStable = joined.all_cases_stable;
isLowSSE = joined.SSE <= trackingCostQ1;
isLowSSdU = joined.SSdU <= inputUseCostQ1;
keepMask = isStable & isLowSSE & isLowSSdU;

stableLowQ = joined(keepMask, :);
if ~isempty(stableLowQ)
    stableLowQ = sortrows(stableLowQ, ["SSE","SSdU"], ["ascend","ascend"]);
end
end


function print_full_run_report(diagnostics, cfg)
caseTable = diagnostics.per_case;
controllerTable = diagnostics.per_controller;
settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));

% Print full controller table (both cases) before summary, with no removals.
allControllersBothCases = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, false);
if isempty(allControllersBothCases)
    fprintf("Controllers table (both case 1 and case 2, no removals): none\n");
else
    fprintf("Controllers table (both case 1 and case 2, no removals):\n");
    allCols = ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag"];
    settleBothCols = allControllersBothCases.Properties.VariableNames(startsWith(allControllersBothCases.Properties.VariableNames, "settle_"));
    allCols = [allCols, settleBothCols, "value_case_stable_c1", "value_case_stable_c2"];
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

if ~isempty(settlingColumns)
    fprintf("Median settling times (h):\n");
    settlingData = table2array(caseTable(:, settlingColumns));
    medianSettling = nanmedian(settlingData, 1);
    for i = 1:numel(settlingColumns)
        fprintf("  %s: %.6g\n", settlingColumns{i}, medianSettling(i));
    end
end

controllerStability = diagnostics.controller_stability;
stableCaseCount = nnz(caseTable.value_case_stable);
stableControllerCount = nnz(controllerStability.all_cases_stable);
fprintf("Possibly stable cases (value nonincreasing): %d/%d\n", stableCaseCount, height(caseTable));
fprintf("Possibly stable controllers (all cases stable): %d/%d\n", stableControllerCount, height(controllerStability));

% Stability counts split by scenario/case id.
case1Mask = (caseTable.case_id == 1);
case2Mask = (caseTable.case_id == 2);
fprintf("Possibly stable case 1 rows: %d/%d\n", nnz(caseTable.value_case_stable & case1Mask), nnz(case1Mask));
fprintf("Possibly stable case 2 rows: %d/%d\n", nnz(caseTable.value_case_stable & case2Mask), nnz(case2Mask));

if ~isempty(diagnostics.missing_pareto_timestamps)
    fprintf("Missing configured Pareto timestamps: %s\n", strjoin(diagnostics.missing_pareto_timestamps, ", "));
end

if isempty(controllerTable)
    return
end

validSSE = controllerTable.SSE(isfinite(controllerTable.SSE));
validSSdU = controllerTable.SSdU(isfinite(controllerTable.SSdU));
if isempty(validSSE) || isempty(validSSdU)
    fprintf("Lower quartiles: unavailable (non-finite controller costs)\n");
    return
end
trackingCostQ1 = prctile(validSSE, 25);
inputUseCostQ1 = prctile(validSSdU, 25);
fprintf("Lower quartiles: SSE <= %.6g, SSdU <= %.6g\n", trackingCostQ1, inputUseCostQ1);

stableControllerList = innerjoin( ...
    controllerStability(controllerStability.all_cases_stable, :), ...
    controllerTable, ...
    "Keys", ["run_label","timestamp"]);
if isempty(stableControllerList)
    fprintf("Stable controllers: none\n");
else
    fprintf("Stable controllers:\n");
    disp(stableControllerList(:, ["run_label","timestamp","SSE","SSdU","J","case_count"]));
end

stableLowQuartileList = diagnostics.stable_low_quartile;
if isempty(stableLowQuartileList)
    fprintf("Stable + lower-quartile (SSE and SSdU): none\n");
else
    fprintf("Stable + lower-quartile (SSE and SSdU):\n");
    disp(stableLowQuartileList(:, ["run_label","timestamp","SSE","SSdU","J","case_count"]));
end

% Controller-level listing: require BOTH cases and no NaN settling times in either case.
if ~isempty(settlingColumns)
    bothCasesTable = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, true);
    if isempty(bothCasesTable)
        fprintf("Controllers table (both cases with complete settling times): none\n");
    else
        settlingBothCols = bothCasesTable.Properties.VariableNames(startsWith(bothCasesTable.Properties.VariableNames, "settle_"));
        fprintf("Controllers table (both case 1 and case 2, excluding any NaN settling time):\n");
        reportCols = ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag", ...
            "value_V1_c1","value_VN_c1","value_V1_c2","value_VN_c2", ...
            settlingBothCols,"value_case_stable_c1","value_case_stable_c2"];
        reportCols = reportCols(ismember(reportCols, bothCasesTable.Properties.VariableNames));
        disp(bothCasesTable(:, reportCols));
    end
end
end


function controllerCaseTable = build_controller_both_case_table(caseTable, controllerTable, settlingColumns, removeNaNSettling)
if isempty(caseTable) || isempty(controllerTable) || isempty(settlingColumns)
    controllerCaseTable = table();
    return
end

case1Table = caseTable(caseTable.case_id == 1, ["run_label","timestamp",settlingColumns,"value_case_stable","value_V1","value_VN"]);
case2Table = caseTable(caseTable.case_id == 2, ["run_label","timestamp",settlingColumns,"value_case_stable","value_V1","value_VN"]);
if isempty(case1Table) || isempty(case2Table)
    controllerCaseTable = table();
    return
end

case1Table = renamevars(case1Table, [settlingColumns, "value_case_stable","value_V1","value_VN"], ...
    [settlingColumns + "_c1", "value_case_stable_c1","value_V1_c1","value_VN_c1"]);
case2Table = renamevars(case2Table, [settlingColumns, "value_case_stable","value_V1","value_VN"], ...
    [settlingColumns + "_c2", "value_case_stable_c2","value_V1_c2","value_VN_c2"]);

controllerCaseTable = innerjoin(case1Table, case2Table, "Keys", ["run_label","timestamp"]);
controllerCaseTable = innerjoin( ...
    controllerCaseTable, ...
    controllerTable(:, ["run_label","timestamp","p","m","SSE","SSdU","J","Q_diag","R1_diag","R2_diag"]), ...
    "Keys", ["run_label","timestamp"]);

if removeNaNSettling
    settlingBothCols = controllerCaseTable.Properties.VariableNames(startsWith(controllerCaseTable.Properties.VariableNames, "settle_"));
    validRowsMask = ~any(isnan(table2array(controllerCaseTable(:, settlingBothCols))), 2);
    controllerCaseTable = controllerCaseTable(validRowsMask, :);
end

if ~isempty(controllerCaseTable)
    controllerCaseTable = sortrows(controllerCaseTable, ["SSE","SSdU"], ["ascend","ascend"]);
end
end


function plot_summary_boxplots(diagnostics)
caseTable = diagnostics.per_case;
controllerTable = diagnostics.per_controller;
if isempty(caseTable)
    return
end

figure("Color","w","Name","test_run_metric_summary");
sgtitle("NMPC Test-Run Diagnostics");

settlingColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "settle_x"));
if ~isempty(settlingColumns)
    subplot(2, 3, 1);
    stateLabels = regexprep(cellstr(settlingColumns), "^settle_(x\\d+)_h$", "$1");
    boxplot(table2array(caseTable(:, settlingColumns)), "Labels", stateLabels);
    xlabel("State");
    ylabel("Settling time (h)");
    title("Settling Time");
    grid off; box off;
end

inputVariationColumns = caseTable.Properties.VariableNames(startsWith(caseTable.Properties.VariableNames, "tot_dU_u"));
if ~isempty(inputVariationColumns)
    subplot(2, 3, 2);
    inputLabels = regexprep(cellstr(inputVariationColumns), "^tot_dU_(u\\d+)$", "$1");
    boxplot(table2array(caseTable(:, inputVariationColumns)), "Labels", inputLabels);
    xlabel("Input");
    ylabel("Total variation (sum |du|)");
    title("Input Variation");
    grid off; box off;
end

subplot(2, 3, 3);
boxplot([caseTable.value_V1, caseTable.value_VN, caseTable.value_noninc_ratio], ...
    "Labels", {'V(1)','V(N)','ratio(DeltaV<=0)'});
xlabel("Metric");
ylabel("Value");
title("Approximate Value Metrics");
grid off; box off;

if ~isempty(controllerTable)
    subplot(2, 3, 4);
    boxplot(controllerTable.SSE, "Labels", {'SSE'});
    set(gca, "YScale", "log");
    ylabel("SSE (log scale)");
    title("Tracking Cost SSE");
    grid off; box off;

    subplot(2, 3, 5);
    boxplot(controllerTable.SSdU, "Labels", {'SSdU'});
    set(gca, "YScale", "log");
    ylabel("SSdU (log scale)");
    title("Input-Use Cost SSdU");
    grid off; box off;

    subplot(2, 3, 6);
    boxplot(controllerTable.J, "Labels", {'J'});
    set(gca, "YScale", "log");
    ylabel("J (log scale)");
    title("Combined Cost J");
    grid off; box off;
end
end
