%% ANALYZE_TEST_RUN_METRICS
% Offline diagnostics for NMPC full-horizon replay data in results/test_run.
% No live simulation is executed.
%
% How to use:
%   1) Ensure replay files exist in:
%        results/test_run/run1_full_f1_no_noise/out_full_*.mat
%        results/test_run/run2_full_f1_no_noise/out_full_*.mat
%   2) Configure:
%        cfg.final_pareto_timestamps  - timestamps to highlight in report
%        cfg.settling_tol             - tolerance (default 2 %)
%        cfg.debug_timestamp          - single timestamp for debug plots
%   3) Run script and inspect:
%        - command window report (aggregate + Pareto-only tables)
%        - debug figure for cfg.debug_timestamp (x1, V, dV per case)
%
% Per-case output columns:
%   - settle_x*_h: settling time per state (hours)
%   - final_err_x*_pct, peak_err_x*_pct: relative tracking errors (%)
%   - tot_dU_u*: total input variation per input (sum |du|)
%   - lyap_V0, lyap_Vf: faux-Lyapunov start/end value
%   - lyap_dV_mean, lyap_dV_max: mean/max increment of faux-Lyapunov
%   - lyap_noninc_ratio: fraction of steps with dV <= 0
%
% Settling-time definition used below:
%   For each state i, relative error is
%       e_rel_i(k) = |y_i(k) - r_i| / max(|r_i|, eps_ref)
%   with r_i the final setpoint value and eps_ref = 1e-9.
%   Settling time t_s,i is the smallest time such that
%       e_rel_i(t) <= 0.02 for all t >= t_s,i.
%   If e_rel_i at the final sample is > 0.02, settling time is NaN.
%
% Additional faux-Lyapunov diagnostic:
%   V(k)  = sum_i (y_i(k) - r_i)^2
%   dV(k) = V(k) - V(k-1)

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
cfg.settling_tol = 0.02;
cfg.debug_timestamp = "20260201_232106";

diagnostics = run_full_run_diagnostics(cfg);
print_full_run_report(diagnostics, cfg);
debug_timestamp_analysis(cfg);

%% FUNCTIONS
function diagnostics = run_full_run_diagnostics(cfg)
%RUN_FULL_RUN_DIAGNOSTICS Load all replay MAT files and build case-wise diagnostics table.
%
% Returns struct:
%   diagnostics.per_case                   table with one row per case
%   diagnostics.mat_file_count             number of MAT files scanned
%   diagnostics.pareto_cases               subset of per_case by timestamp
%   diagnostics.missing_pareto_timestamps  configured timestamps not found
arguments
    cfg struct
end

rows = [];
matCount = 0;

for k = 1:numel(cfg.run_folders)
    runLabel = cfg.run_folders(k).label;
    runDir = fullfile(cfg.test_run_root, cfg.run_folders(k).subfolder);
    if ~isfolder(runDir)
        warning("Run folder missing: %s", runDir);
        continue
    end

    matList = dir(fullfile(runDir, "out_full_*.mat"));
    matCount = matCount + numel(matList);

    for m = 1:numel(matList)
        matPath = fullfile(runDir, matList(m).name);
        S = load(matPath, "ts", "out");
        if ~isfield(S, "ts") || ~isfield(S, "out") || ~isfield(S.out, "case")
            warning("Skipping %s (missing ts/out/case).", matPath);
            continue
        end

        ts = string(S.ts);
        out = S.out;

        for c = 1:numel(out.case)
            ci = out.case(c);
            metrics = summarize_case_metrics(out, ci, cfg.settling_tol);

            row = struct();
            row.run_label = string(runLabel);
            row.timestamp = ts;
            row.case_id = ci.case_id;
            if isfield(out, "tf")
                row.tf_h = double(out.tf);
            else
                row.tf_h = NaN;
            end
            row.N_steps = size(ci.Y, 1);

            for s = 1:numel(metrics.settling_times_h)
                row.(sprintf("settle_x%d_h", s)) = metrics.settling_times_h(s);
                row.(sprintf("final_err_x%d_pct", s)) = metrics.final_error_pct(s);
                row.(sprintf("peak_err_x%d_pct", s)) = metrics.peak_error_pct(s);
            end

            for u = 1:numel(metrics.total_input_variation)
                row.(sprintf("tot_dU_u%d", u)) = metrics.total_input_variation(u);
            end

            row.lyap_V0 = metrics.lyap_V0;
            row.lyap_Vf = metrics.lyap_Vf;
            row.lyap_dV_mean = metrics.lyap_dV_mean;
            row.lyap_dV_max = metrics.lyap_dV_max;
            row.lyap_noninc_ratio = metrics.lyap_noninc_ratio;

            rows = [rows; row]; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    perCase = table();
else
    perCase = struct2table(rows);
end

diagnostics = struct();
diagnostics.per_case = perCase;
diagnostics.mat_file_count = matCount;

if isempty(perCase)
    diagnostics.pareto_cases = table();
    diagnostics.missing_pareto_timestamps = cfg.final_pareto_timestamps;
    return
end

diagnostics.pareto_cases = perCase(ismember(perCase.timestamp, cfg.final_pareto_timestamps), :);
foundTs = unique(perCase.timestamp);
diagnostics.missing_pareto_timestamps = setdiff(cfg.final_pareto_timestamps, foundTs);
end


function metrics = summarize_case_metrics(outStruct, caseStruct, settlingTol)
%SUMMARIZE_CASE_METRICS Compute all per-case scalar diagnostics.
arguments
    outStruct struct
    caseStruct struct
    settlingTol (1,1) double {mustBePositive} = 0.02
end

Y = double(caseStruct.Y);
Ysp = double(caseStruct.Ysp);
U = double(caseStruct.U);

if isfield(caseStruct, "dt")
    dt_h = double(caseStruct.dt);
elseif isfield(outStruct, "dt")
    dt_h = double(outStruct.dt);
else
    error("Missing dt in case/out struct.");
end
t_h = (0:size(Y,1)-1).' * dt_h;

[settlingTimes_h, finalErrPct, peakErrPct] = compute_settling_times(Y, Ysp, t_h, settlingTol);
totalVar = compute_total_input_variation(U);
[V, dV] = compute_faux_lyapunov(Y, Ysp);

metrics = struct();
metrics.settling_times_h = settlingTimes_h;
metrics.final_error_pct = finalErrPct;
metrics.peak_error_pct = peakErrPct;
metrics.total_input_variation = totalVar;
metrics.lyap_V0 = V(1);
metrics.lyap_Vf = V(end);
metrics.lyap_dV_mean = mean(dV(2:end), "omitnan");
metrics.lyap_dV_max = max(dV(2:end), [], "omitnan");
metrics.lyap_noninc_ratio = mean(dV(2:end) <= 0, "omitnan");
end


function [settlingTimes_h, finalErrPct, peakErrPct] = compute_settling_times(Y, Ysp, t_h, tol)
%COMPUTE_SETTLING_TIMES Compute settling time in hours using "stay within band forever" rule.
%
% Notes:
%   - Relative error normalization uses |reference| (with eps floor).
%   - If final relative error is outside tolerance, settling time is NaN.
arguments
    Y double
    Ysp double
    t_h double
    tol (1,1) double {mustBeGreaterThan(tol,0)} = 0.02
end

nStates = size(Y, 2);
settlingTimes_h = NaN(1, nStates);
finalErrPct = NaN(1, nStates);
peakErrPct = NaN(1, nStates);
eps_ref = 1e-9;

for s = 1:nStates
    refVal = Ysp(end, s);
    scale = max(abs(refVal), eps_ref);

    errRel = abs(Y(:, s) - refVal) ./ scale;
    finalErrPct(s) = 100 * errRel(end);
    peakErrPct(s) = 100 * max(errRel);

    if finalErrPct(s) > 100 * tol
        settlingTimes_h(s) = NaN;
        continue
    end

    % First index from which ALL future samples stay inside the tolerance band.
    futureMax = flipud(cummax(flipud(errRel)));
    idx = find(futureMax <= tol, 1, "first");
    if ~isempty(idx)
        settlingTimes_h(s) = t_h(idx);
    end
end
end


function totalVar = compute_total_input_variation(U)
%COMPUTE_TOTAL_INPUT_VARIATION Return per-input total variation sum(|du|).
if isempty(U)
    totalVar = [];
    return
end
if size(U,1) < 2
    totalVar = zeros(1, size(U,2));
    return
end
dU = diff(U, 1, 1);
totalVar = sum(abs(dU), 1);
end


function [V, dV] = compute_faux_lyapunov(Y, Ysp)
%COMPUTE_FAUX_LYAPUNOV Compute V(k)=||y-r||^2 and dV(k)=V(k)-V(k-1).
E = Y - Ysp;
V = sum(E.^2, 2);
dV = [NaN; diff(V)];
end


function print_full_run_report(diagnostics, cfg)
%PRINT_FULL_RUN_REPORT Print concise aggregate and Pareto-focused diagnostics.
perCase = diagnostics.per_case;

fprintf("=== NMPC Full-Horizon Replay Report ===\n");
fprintf("Run folders scanned: %s\n", strjoin(string({cfg.run_folders.subfolder}), ", "));
fprintf("MAT-files detected: %d\n", diagnostics.mat_file_count);
fprintf("Case rows extracted: %d\n", height(perCase));

if isempty(perCase)
    fprintf("No case data found.\n");
    return
end

settleCols = perCase.Properties.VariableNames(startsWith(perCase.Properties.VariableNames, "settle_x"));
if ~isempty(settleCols)
    medSettle = varfun(@nanmedian, perCase(:, settleCols));
    fprintf("\nMedian settling times across all cases (hours):\n");
    disp(medSettle);
end

fprintf("Pareto timestamps requested: %d\n", numel(cfg.final_pareto_timestamps));
paretoCases = diagnostics.pareto_cases;
if isempty(paretoCases)
    fprintf("Pareto timestamps present: 0\n");
else
    fprintf("Pareto timestamps present: %d\n", numel(unique(paretoCases.timestamp)));
end
if ~isempty(diagnostics.missing_pareto_timestamps)
    fprintf("Missing Pareto controllers: %s\n", strjoin(diagnostics.missing_pareto_timestamps, ", "));
end

if isempty(paretoCases)
    return
end

keyCols = ["run_label","timestamp","case_id","tf_h","N_steps"];
stateCols = perCase.Properties.VariableNames(startsWith(perCase.Properties.VariableNames, "settle_x"));
errCols = perCase.Properties.VariableNames(contains(perCase.Properties.VariableNames, "_err_"));
uCols = perCase.Properties.VariableNames(startsWith(perCase.Properties.VariableNames, "tot_dU_u"));
lyapCols = ["lyap_V0","lyap_Vf","lyap_dV_mean","lyap_dV_max","lyap_noninc_ratio"];
showCols = [keyCols, stateCols, errCols, uCols, lyapCols];
showCols = showCols(ismember(showCols, paretoCases.Properties.VariableNames));

fprintf("\nPareto controller diagnostics (per case):\n");
disp(paretoCases(:, showCols));

if ~isempty(uCols)
    fprintf("Median total input variation over Pareto cases:\n");
    disp(varfun(@nanmedian, paretoCases(:, uCols)));
end

fprintf("Median faux-Lyapunov stats over Pareto cases:\n");
disp(varfun(@nanmedian, paretoCases(:, lyapCols)));
end


function debug_timestamp_analysis(cfg)
%DEBUG_TIMESTAMP_ANALYSIS Plot and print detailed diagnostics for one timestamp.
%
% Plots per case:
%   - x1 and x1 setpoint over time
%   - faux-Lyapunov V over time
%   - increment dV over time
if ~isfield(cfg, "debug_timestamp") || strlength(cfg.debug_timestamp) == 0
    return
end

ts = string(cfg.debug_timestamp);
[matPath, runLabel] = find_test_run_mat(ts, cfg);
if strlength(matPath) == 0
    fprintf("\nDebug timestamp %s not found in test_run folders.\n", ts);
    return
end

S = load(matPath, "out");
if ~isfield(S, "out") || ~isfield(S.out, "case")
    fprintf("\nDebug file %s missing out/case.\n", matPath);
    return
end
out = S.out;

fprintf("\n=== Debug timestamp %s (%s) ===\n", ts, runLabel);
fig = figure("Color","w","Name","debug_"+ts);

for c = 1:numel(out.case)
    ci = out.case(c);
    metrics = summarize_case_metrics(out, ci, cfg.settling_tol);
    [V, dV] = compute_faux_lyapunov(double(ci.Y), double(ci.Ysp));

    if isfield(ci, "dt")
        dt_h = double(ci.dt);
    elseif isfield(out, "dt")
        dt_h = double(out.dt);
    else
        dt_h = NaN;
    end
    t_h = (0:size(ci.Y,1)-1).' * dt_h;

    subplot(numel(out.case), 3, (c-1)*3 + 1);
    plot(t_h, ci.Y(:,1), "LineWidth", 1.4); hold on;
    plot(t_h, ci.Ysp(:,1), "--", "LineWidth", 1.2);
    xlabel("t (h)"); ylabel("x1");
    title(sprintf("case %d: x1", ci.case_id));
    legend("x1","x1 sp","Location","best");
    grid off; box off;

    subplot(numel(out.case), 3, (c-1)*3 + 2);
    plot(t_h, V, "LineWidth", 1.4);
    xlabel("t (h)"); ylabel("V");
    title("faux-Lyapunov V");
    grid off; box off;

    subplot(numel(out.case), 3, (c-1)*3 + 3);
    plot(t_h, dV, "LineWidth", 1.2);
    xlabel("t (h)"); ylabel("dV");
    title("V(k)-V(k-1)");
    grid off; box off;

    fprintf("Case %d:\n", ci.case_id);
    fprintf("  settle x1 (h): %.6g\n", metrics.settling_times_h(1));
    fprintf("  final/peak x1 err (%%): %.4g / %.4g\n", metrics.final_error_pct(1), metrics.peak_error_pct(1));
    fprintf("  total dU: %s\n", mat2str(metrics.total_input_variation, 5));
    fprintf("  V0=%.6g, Vf=%.6g, mean(dV)=%.6g, max(dV)=%.6g, noninc_ratio=%.4f\n", ...
        metrics.lyap_V0, metrics.lyap_Vf, metrics.lyap_dV_mean, metrics.lyap_dV_max, metrics.lyap_noninc_ratio);
end
end


function [matPath, runLabel] = find_test_run_mat(ts, cfg)
%FIND_TEST_RUN_MAT Resolve timestamp to replay MAT path and run label.
matPath = "";
runLabel = "";

for k = 1:numel(cfg.run_folders)
    runLabelCandidate = cfg.run_folders(k).label;
    runDir = fullfile(cfg.test_run_root, cfg.run_folders(k).subfolder);
    candidate = fullfile(runDir, "out_full_" + ts + ".mat");
    if isfile(candidate)
        matPath = candidate;
        runLabel = runLabelCandidate;
        return
    end
end
end
