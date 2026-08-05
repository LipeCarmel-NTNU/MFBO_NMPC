%% Audit selected schedule controller against schedule metrics and main Pareto data
% This script is intentionally read-only. It checks whether the selected
% controller comment block is supported by schedule metrics or by the
% original run1/run2 Pareto analysis.

clear; clc

scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
resultsRoot = fullfile(projectRoot, "results");

selectedId = "ts_20260211_122653_modified";
selectedTs = "20260211_122653";

scheduleCsv = fullfile(resultsRoot, "numerical results", "setpoint_schedule_metrics_same_noise.csv");
rawScheduleCsv = fullfile(resultsRoot, "setpoint_schedule_xsp_7_13_16", "same_noise", "results_schedule.csv");
finalParetoCsv = fullfile(resultsRoot, "numerical results", "final_pareto_frontier_f1_noisy_noiseless_metrics.csv");

fprintf("\n=== Selected schedule pick audit ===\n");
fprintf("selected controller_id: %s\n", selectedId);
fprintf("base timestamp:         %s\n", selectedTs);

Tsch = readtable(scheduleCsv, "TextType", "string");
Traw = readtable(rawScheduleCsv, "TextType", "string");

sel = require_one(Tsch, string(Tsch.controller_id) == selectedId, "selected schedule row");
bench = require_one(Tsch, logical(Tsch.is_benchmark), "schedule benchmark row");

fprintf("\n--- Schedule metrics source ---\n");
fprintf("schedule CSV: %s\n", scheduleCsv);
fprintf("selected row source: %s\n", string(sel.source));
fprintf("selected row timestamp: %s\n", string(sel.timestamp));
fprintf("benchmark row: %s\n", string(bench.controller_id));

fprintf("\n--- Selected vs benchmark, schedule xsp_7_13_16 same_noise ---\n");
print_ratio("SSE_total", double(sel.SSE_total), double(bench.SSE_total), false);
print_ratio("SSdU_total", double(sel.SSdU_total), double(bench.SSdU_total), false);
print_ratio("wall_time_s", double(sel.wall_time_s), double(bench.wall_time_s), false);
for stateIdx = 1:3
    selIae = state_iae_total(sel, stateIdx);
    benchIae = state_iae_total(bench, stateIdx);
    print_ratio(sprintf("x%d_IAE_total", stateIdx), selIae, benchIae, false);
end

fprintf("\nRaw schedule CSV contains selected/unmodified duplicates before analysis drop:\n");
rawMatch = Traw(contains(string(Traw.controller_id), "ts_20260211_122653"), ...
    intersect(["controller_id","source","timestamp","SSE","SSdU","runtime_s"], string(Traw.Properties.VariableNames), "stable"));
disp(rawMatch);

fprintf("\n--- Main Pareto reconstruction from results/run1 and results/run2 ---\n");
TallOpt = table();
runKeys = ["run1", "run2"];
for k = 1:numel(runKeys)
    runKey = runKeys(k);
    Trun = readtable(fullfile(resultsRoot, runKey, "results.csv"), "TextType", "string");
    Trun.run_key = repmat(runKey, height(Trun), 1);
    Trun.iteration = (1:height(Trun)).';
    Trun.is_optimization = double(Trun.iteration) > 20;
    Trun.is_run_pareto = false(height(Trun), 1);

    optMask = logical(Trun.is_optimization);
    runParetoMask = compute_pareto_mask_local(double(Trun.SSE(optMask)), double(Trun.SSdU(optMask)));
    optIdx = find(optMask);
    Trun.is_run_pareto(optIdx(runParetoMask)) = true;

    TallOpt = [TallOpt; Trun(optMask, :)]; %#ok<AGROW>

    targetRunMask = string(Trun.timestamp) == selectedTs;
    if any(targetRunMask)
        targetRow = Trun(targetRunMask, :);
        fprintf("%s target: iteration=%d, optimization=%d, run-pareto=%d, SSE=%.15g, SSdU=%.15g, f=%.15g\n", ...
            runKey, double(targetRow.iteration), double(targetRow.is_optimization), double(targetRow.is_run_pareto), ...
            double(targetRow.SSE), double(targetRow.SSdU), double(targetRow.theta_1));
        report_ranks(Trun(optMask, :), targetRow, "within " + runKey + " optimization rows");
        report_ranks(Trun(logical(Trun.is_run_pareto), :), targetRow, "within " + runKey + " run Pareto rows");
    end
end

combinedParetoMask = compute_pareto_mask_local(double(TallOpt.SSE), double(TallOpt.SSdU));
TallOpt.is_combined_pareto = combinedParetoMask;
targetCombined = require_one(TallOpt, string(TallOpt.timestamp) == selectedTs & string(TallOpt.run_key) == "run2", "combined target row");
fprintf("combined target: combined-pareto=%d, SSE=%.15g, SSdU=%.15g\n", ...
    double(targetCombined.is_combined_pareto), double(targetCombined.SSE), double(targetCombined.SSdU));
report_ranks(TallOpt, targetCombined, "within combined optimization rows");
report_ranks(TallOpt(combinedParetoMask, :), targetCombined, "within combined main Pareto rows");

fprintf("\nLowest-SSE combined main Pareto rows:\n");
showCols = intersect(["run_key","timestamp","SSE","SSdU","J","theta_1","iteration"], string(TallOpt.Properties.VariableNames), "stable");
Tfront = sortrows(TallOpt(combinedParetoMask, showCols), "SSE", "ascend");
disp(Tfront(1:min(8, height(Tfront)), :));

if isfile(finalParetoCsv)
    fprintf("\n--- Final f=1 fidelity table written by main_pareto ---\n");
    Tf1 = readtable(finalParetoCsv, "TextType", "string");
    targetF1 = require_one(Tf1, string(Tf1.timestamp) == selectedTs & string(Tf1.run_key) == "run2", "final f1 target row");
    report_metric_rank(Tf1, targetF1, "noisy_Jtrack", "within final f=1 Pareto table");
    report_metric_rank(Tf1, targetF1, "noisy_JTV", "within final f=1 Pareto table");
    fprintf("target final f=1 noisy_Jtrack=%.15g, noisy_JTV=%.15g\n", ...
        double(targetF1.noisy_Jtrack), double(targetF1.noisy_JTV));
end

function row = require_one(T, mask, label)
idx = find(mask);
if isempty(idx)
    error("Could not find %s.", label);
end
if numel(idx) > 1
    warning("Found %d rows for %s; using the first.", numel(idx), label);
end
row = T(idx(1), :);
end

function print_ratio(label, selectedValue, benchmarkValue, higherIsBetter)
ratio = selectedValue / benchmarkValue;
pct = 100 * ratio;
if higherIsBetter
    deltaText = sprintf("%.4f%% higher", 100 * (ratio - 1));
else
    if ratio <= 1
        deltaText = sprintf("%.4f%% lower", 100 * (1 - ratio));
    else
        deltaText = sprintf("%.4f%% higher", 100 * (ratio - 1));
    end
end
fprintf("%-14s selected=%.15g benchmark=%.15g ratio=%.6f (%.4f%% of benchmark, %s)\n", ...
    label, selectedValue, benchmarkValue, ratio, pct, deltaText);
end

function v = state_iae_total(T, stateIdx)
v = double(T.(sprintf("IAE_x%d_c1", stateIdx))) + double(T.(sprintf("IAE_x%d_c2", stateIdx)));
end

function isPareto = compute_pareto_mask_local(Jtrack, Jtv)
n = numel(Jtrack);
isPareto = true(n, 1);
for i = 1:n
    dominated = (Jtrack <= Jtrack(i)) & (Jtv <= Jtv(i)) & ...
        ((Jtrack < Jtrack(i)) | (Jtv < Jtv(i)));
    dominated(i) = false;
    if any(dominated)
        isPareto(i) = false;
    end
end
end

function report_ranks(T, targetRow, scopeLabel)
if isempty(T)
    fprintf("rank scope empty: %s\n", scopeLabel);
    return
end
report_metric_rank(T, targetRow, "SSE", scopeLabel);
report_metric_rank(T, targetRow, "SSdU", scopeLabel);
end

function report_metric_rank(T, targetRow, metricName, scopeLabel)
values = double(T.(metricName));
targetValue = double(targetRow.(metricName));
[~, ord] = sort(values, "ascend");
rank = find(abs(values(ord) - targetValue) <= max(1e-12, 1e-12 * abs(targetValue)), 1, "first");
fprintf("%s rank by lowest %s: %d of %d\n", scopeLabel, metricName, rank, height(T));
end
