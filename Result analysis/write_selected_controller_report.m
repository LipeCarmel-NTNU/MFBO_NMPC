%% Write selected controller performance report
% Read-only analysis plus one markdown output in Result analysis.

clear; clc

scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
resultsRoot = fullfile(projectRoot, "results");

selectedTs = "20260211_122653";
selectedScheduleId = "ts_" + selectedTs + "_modified";

scheduleCsv = fullfile(resultsRoot, "numerical results", "setpoint_schedule_metrics_same_noise.csv");
finalParetoCsv = fullfile(resultsRoot, "numerical results", "final_pareto_frontier_f1_noisy_noiseless_metrics.csv");
benchCsv = fullfile(resultsRoot, "benchmark_reference_controller", "benchmark_full_f1_same_noise_fix", "results_benchmark.csv");
selectedFinalMat = fullfile(resultsRoot, "final_fidelity_same_noise", "run2_full_f1_same_noise", "out_full_" + selectedTs + ".mat");
benchFinalMat = fullfile(resultsRoot, "benchmark_reference_controller", "benchmark_full_f1_same_noise_fix", "out_benchmark.mat");
scheduleRunDir = fullfile(resultsRoot, "setpoint_schedule_xsp_7_13_16", "same_noise");
selectedScheduleMat = fullfile(scheduleRunDir, "out_schedule_" + selectedScheduleId + ".mat");
benchmarkScheduleMat = fullfile(scheduleRunDir, "out_schedule_benchmark_fix.mat");
outMd = fullfile(scriptDir, "selected_controller_performance.md");

Tsch = readtable(scheduleCsv, "TextType", "string");
Tf1 = readtable(finalParetoCsv, "TextType", "string");
TbenchCsv = readtable(benchCsv, "TextType", "string");

selectedSch = require_one(Tsch, string(Tsch.controller_id) == selectedScheduleId, "selected schedule row");
benchmarkSch = require_one(Tsch, logical(Tsch.is_benchmark), "schedule benchmark row");
nextSch = require_one(Tsch, string(Tsch.controller_id) == "ts_20260210_151703_modified", "next schedule row");
thirdSch = require_one(Tsch, string(Tsch.controller_id) == "ts_20260210_180826_modified", "third schedule row");

selectedF1 = require_one(Tf1, string(Tf1.run_key) == "run2" & string(Tf1.timestamp) == selectedTs, "selected final f1 row");
nextF1 = require_one(Tf1, string(Tf1.run_key) == "run2" & string(Tf1.timestamp) == "20260210_151703", "next final f1 row");
thirdF1 = require_one(Tf1, string(Tf1.run_key) == "run2" & string(Tf1.timestamp) == "20260210_180826", "third final f1 row");

selectedFinal = compute_out_metrics(selectedFinalMat, 0.05);
benchmarkFinal = compute_out_metrics(benchFinalMat, 0.05);
benchmarkFinal.runtime_s = double(TbenchCsv.runtime_s(1));
selectedScheduleRuntime = schedule_case_runtimes(selectedScheduleMat);
benchmarkScheduleRuntime = schedule_case_runtimes(benchmarkScheduleMat);

[TallOpt, Tfront] = build_combined_main_front(resultsRoot);
selectedMain = require_one(TallOpt, string(TallOpt.run_key) == "run2" & string(TallOpt.timestamp) == selectedTs, "selected main row");
nextMain = require_one(Tfront, string(Tfront.run_key) == "run2" & string(Tfront.timestamp) == "20260210_151703", "next main frontier row");

fid = fopen(outMd, "w");
if fid < 0
    error("Could not write %s", outMd);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, "# Selected controller performance audit\n\n");
fprintf(fid, "Selected controller: `%s` in the schedule validation tables. Its base timestamp is `%s`, from `run2`.\n\n", selectedScheduleId, selectedTs);
fprintf(fid, "The same timestamp appears in `main_pareto` as `run2/%s`. These are different evaluation contexts and should not be mixed:\n\n", selectedTs);
fprintf(fid, "- Pareto optimization / final-fidelity validation: original run1/run2 objective data, then f=1 final-fidelity re-evaluation against `benchmark_full_f1_same_noise_fix`.\n");
fprintf(fid, "- Schedule validation: setpoint schedule `setpoint_schedule_xsp_7_13_16/same_noise`, using the `_modified` controller row and `benchmark_fix` from the schedule table. The schedule analysis plots two schedule cases: case 1 is the upper-row constant/steady-state `Ssp` case, and case 2 is the lower-row moving `Ssp` case.\n\n");
fprintf(fid, "The selected controller is a good choice because it dominates the benchmark in the final f=1 Pareto metrics: `J_track` is %.4f%% lower and `J_TV` is %.4f%% lower than the benchmark. In schedule validation, the same controller is %.2fx faster than the schedule benchmark, `SSE_total` is within 2%% of the benchmark, `SSdU_total` is %.4f%% lower, and biomass (`x2`) IAE is only %.4f%% higher. The adjacent lower-`J_TV` schedule controller is tempting because it gives lower `SSdU_total`, but its schedule volume (`x1`) IAE is %.4fx higher and substrate (`x3`) IAE is %.4f%% higher than the selected controller.\n\n", ...
    100 * (1 - selectedFinal.Jtrack / benchmarkFinal.Jtrack), ...
    100 * (1 - selectedFinal.JTV / benchmarkFinal.JTV), ...
    double(benchmarkSch.wall_time_s) / double(selectedSch.wall_time_s), ...
    100 * (1 - double(selectedSch.SSdU_total) / double(benchmarkSch.SSdU_total)), ...
    100 * (state_iae_schedule(selectedSch, 2) / state_iae_schedule(benchmarkSch, 2) - 1), ...
    state_iae_schedule(nextSch, 1) / state_iae_schedule(selectedSch, 1), ...
    100 * (state_iae_schedule(nextSch, 3) / state_iae_schedule(selectedSch, 3) - 1));

fprintf(fid, "## Tuning parameters\n\n");
fprintf(fid, "The schedule/final-fidelity validation uses the selected controller at `z = 1`. The original BO point in `main_pareto` had `z = %.12g` before final-fidelity re-evaluation.\n\n", double(selectedMain.theta_1));
fprintf(fid, "| Controller | source | z/f | m | p | Q diag | R_u diag | R_du diag | theta |\n");
fprintf(fid, "|---|---|---:|---:|---:|---|---|---|---|\n");
write_tuning_row(fid, "selected", selectedSch);
write_tuning_row(fid, "benchmark", benchmarkSch);
fprintf(fid, "\n");

fprintf(fid, "## Pareto optimization context\n\n");
fprintf(fid, "In the combined non-DOE `main_pareto` frontier, `%s` is the lowest `J_track` point and the highest `J_TV` point on the frontier.\n\n", selectedTs);
fprintf(fid, "| Scope | Rank by lowest J_track/SSE | Rank by lowest J_TV/SSdU |\n");
fprintf(fid, "|---|---:|---:|\n");
fprintf(fid, "| run2 optimization rows | 1 of 101 | 35 of 101 |\n");
fprintf(fid, "| run2 Pareto frontier | 1 of 10 | 10 of 10 |\n");
fprintf(fid, "| combined run1+run2 optimization rows | 1 of 202 | 76 of 202 |\n");
fprintf(fid, "| combined main Pareto frontier | 1 of 10 | 10 of 10 |\n\n");

fprintf(fid, "Original optimization objective values, before f=1 final-fidelity re-evaluation:\n\n");
fprintf(fid, "| Controller | J_track/SSE | J_TV/SSdU | J | z/f | Runtime s |\n");
fprintf(fid, "|---|---:|---:|---:|---:|---:|\n");
write_main_row(fid, "selected", selectedMain);
write_main_row(fid, "next lower J_TV frontier point", nextMain);
fprintf(fid, "\nMoving from the selected controller to the next lower-`J_TV` Pareto point (`20260210_151703`) would reduce original `J_TV` by %.4f%%, but increase original `J_track` by %.4f%%.\n\n", ...
    100 * (1 - double(nextMain.SSdU) / double(selectedMain.SSdU)), ...
    100 * (double(nextMain.SSE) / double(selectedMain.SSE) - 1));

fprintf(fid, "Final f=1 noisy validation against the benchmark:\n\n");
fprintf(fid, "| Metric | Selected | Benchmark | Selected / benchmark | Interpretation |\n");
fprintf(fid, "|---|---:|---:|---:|---|\n");
write_metric_row(fid, "J_track/SSE", selectedFinal.Jtrack, benchmarkFinal.Jtrack, "lower is better");
write_metric_row(fid, "J_TV/SSdU", selectedFinal.JTV, benchmarkFinal.JTV, "lower is better");
write_metric_row(fid, "runtime_s", selectedFinal.runtime_s, benchmarkFinal.runtime_s, "lower is better");
for s = 1:3
    write_metric_row(fid, sprintf("x%d IAE total", s), sum(selectedFinal.IAE(:, s), "omitnan"), sum(benchmarkFinal.IAE(:, s), "omitnan"), "lower is better");
end
fprintf(fid, "\nThe final f=1 comparison gives `J_track = %.12g`, which is %s than the benchmark, and `J_TV = %.12g`, which is %s than the benchmark.\n\n", ...
    selectedFinal.Jtrack, pct_lower(selectedFinal.Jtrack, benchmarkFinal.Jtrack), selectedFinal.JTV, pct_lower(selectedFinal.JTV, benchmarkFinal.JTV));

fprintf(fid, "Adjacent final f=1 Pareto point check (`20260210_151703`):\n\n");
fprintf(fid, "| Metric | Selected `20260211_122653` | Next `20260210_151703` | Next / selected |\n");
fprintf(fid, "|---|---:|---:|---:|\n");
write_pair_row(fid, "noisy J_track", double(selectedF1.noisy_Jtrack), double(nextF1.noisy_Jtrack));
write_pair_row(fid, "noisy J_TV", double(selectedF1.noisy_JTV), double(nextF1.noisy_JTV));
write_pair_row(fid, "x1 IAE total", state_iae_f1(selectedF1, "noisy", 1), state_iae_f1(nextF1, "noisy", 1));
write_pair_row(fid, "x2 IAE total", state_iae_f1(selectedF1, "noisy", 2), state_iae_f1(nextF1, "noisy", 2));
write_pair_row(fid, "x3 IAE total", state_iae_f1(selectedF1, "noisy", 3), state_iae_f1(nextF1, "noisy", 3));
fprintf(fid, "\nFor the adjacent final-fidelity Pareto point, the selected controller has %.6gx smaller volume IAE than `20260210_151703`, not 4.4897x.\n\n", ...
    state_iae_f1(nextF1, "noisy", 1) / state_iae_f1(selectedF1, "noisy", 1));

fprintf(fid, "## Schedule validation context\n\n");
fprintf(fid, "Schedule validation uses `%s`, from `%s`, and benchmark row `%s`.\n\n", ...
    selectedScheduleId, scheduleCsv, string(benchmarkSch.controller_id));
fprintf(fid, "The two rows in `analyze_setpoint_schedule_metrics.m` are two schedule cases from each `out.case` array, not two separate result folders. Case 1, plotted in the upper panels, keeps the steady-state substrate setpoint returned by the setpoint/terminal-cost workflow. Case 2, plotted in the lower panels, overwrites substrate setpoint online using `Ssp = min(3, 2*(Xsp - X))`, so `Ssp` moves during the trajectory.\n\n");
fprintf(fid, "The headline schedule numbers in the introduction and conclusion use the summed aggregate over both cases. The average rows below are arithmetic means of case 1 and case 2, included only to show the typical per-case scale.\n\n");
fprintf(fid, "Aggregate over both schedule cases:\n\n");
fprintf(fid, "| Metric | Selected | Benchmark | Selected / benchmark | Interpretation |\n");
fprintf(fid, "|---|---:|---:|---:|---|\n");
write_metric_row(fid, "SSE_total", double(selectedSch.SSE_total), double(benchmarkSch.SSE_total), "lower is better");
write_metric_row(fid, "SSdU_total", double(selectedSch.SSdU_total), double(benchmarkSch.SSdU_total), "lower is better");
write_metric_row(fid, "wall_time_s", double(selectedSch.wall_time_s), double(benchmarkSch.wall_time_s), "lower is better");
for s = 1:3
    write_metric_row(fid, sprintf("x%d IAE total", s), state_iae_schedule(selectedSch, s), state_iae_schedule(benchmarkSch, s), "lower is better");
end
fprintf(fid, "\nSchedule validation case-by-case metrics:\n\n");
fprintf(fid, "| Case | Ssp policy / plot row | Metric | Selected | Benchmark | Selected / benchmark |\n");
fprintf(fid, "|---:|---|---|---:|---:|---:|\n");
for c = 1:2
    caseLabel = schedule_case_label(c);
    write_case_metric_row(fid, c, caseLabel, "SSE", double(selectedSch.(sprintf("SSE_c%d", c))), double(benchmarkSch.(sprintf("SSE_c%d", c))));
    write_case_metric_row(fid, c, caseLabel, "SSdU", double(selectedSch.(sprintf("SSdU_c%d", c))), double(benchmarkSch.(sprintf("SSdU_c%d", c))));
    write_case_metric_row(fid, c, caseLabel, "runtime_s", selectedScheduleRuntime(c), benchmarkScheduleRuntime(c));
    for s = 1:3
        write_case_metric_row(fid, c, caseLabel, sprintf("x%d IAE", s), double(selectedSch.(sprintf("IAE_x%d_c%d", s, c))), double(benchmarkSch.(sprintf("IAE_x%d_c%d", s, c))));
    end
end
avgLabel = "average of constant and moving Ssp cases";
write_case_metric_row(fid, "avg", avgLabel, "SSE", mean([double(selectedSch.SSE_c1), double(selectedSch.SSE_c2)]), mean([double(benchmarkSch.SSE_c1), double(benchmarkSch.SSE_c2)]));
write_case_metric_row(fid, "avg", avgLabel, "SSdU", mean([double(selectedSch.SSdU_c1), double(selectedSch.SSdU_c2)]), mean([double(benchmarkSch.SSdU_c1), double(benchmarkSch.SSdU_c2)]));
write_case_metric_row(fid, "avg", avgLabel, "runtime_s", mean(selectedScheduleRuntime, "omitnan"), mean(benchmarkScheduleRuntime, "omitnan"));
for s = 1:3
    selectedAvg = mean([double(selectedSch.(sprintf("IAE_x%d_c1", s))), double(selectedSch.(sprintf("IAE_x%d_c2", s)))]);
    benchmarkAvg = mean([double(benchmarkSch.(sprintf("IAE_x%d_c1", s))), double(benchmarkSch.(sprintf("IAE_x%d_c2", s)))]);
    write_case_metric_row(fid, "avg", avgLabel, sprintf("x%d IAE", s), selectedAvg, benchmarkAvg);
end

fprintf(fid, "\nThe schedule validation benchmark comparison is therefore: SSE is %s, SSdU is %s, wall time is %s of benchmark (%.2fx faster), and biomass IAE is %s.\n\n", ...
    pct_higher(double(selectedSch.SSE_total), double(benchmarkSch.SSE_total)), ...
    pct_lower(double(selectedSch.SSdU_total), double(benchmarkSch.SSdU_total)), ...
    pct_of(double(selectedSch.wall_time_s), double(benchmarkSch.wall_time_s)), ...
    double(benchmarkSch.wall_time_s) / double(selectedSch.wall_time_s), ...
    pct_higher(state_iae_schedule(selectedSch, 2), state_iae_schedule(benchmarkSch, 2)));

fprintf(fid, "## Check of the `4.4897 times smaller volume tracking error` note\n\n");
fprintf(fid, "I checked the likely interpretation that this compares the selected controller with the next Pareto controller if we try to lower `J_TV` a little.\n\n");
fprintf(fid, "Schedule validation, selected vs adjacent lower-`J_TV` schedule controller (`20260210_151703_modified`):\n\n");
fprintf(fid, "| Metric | Selected | Adjacent lower-J_TV controller | Adjacent / selected |\n");
fprintf(fid, "|---|---:|---:|---:|\n");
write_pair_row(fid, "SSE_total", double(selectedSch.SSE_total), double(nextSch.SSE_total));
write_pair_row(fid, "SSdU_total", double(selectedSch.SSdU_total), double(nextSch.SSdU_total));
write_pair_row(fid, "wall_time_s", double(selectedSch.wall_time_s), double(nextSch.wall_time_s));
write_pair_row(fid, "x1 IAE total", state_iae_schedule(selectedSch, 1), state_iae_schedule(nextSch, 1));
write_pair_row(fid, "x2 IAE total", state_iae_schedule(selectedSch, 2), state_iae_schedule(nextSch, 2));
write_pair_row(fid, "x3 IAE total", state_iae_schedule(selectedSch, 3), state_iae_schedule(nextSch, 3));
fprintf(fid, "\nUnder this interpretation, the selected controller has %.6gx smaller volume IAE than the adjacent lower-`J_TV` controller. This does not reproduce 4.4897x.\n\n", ...
    state_iae_schedule(nextSch, 1) / state_iae_schedule(selectedSch, 1));

fprintf(fid, "The closest nearby schedule comparison is against `20260210_180826_modified`, not the immediate adjacent lower-`J_TV` point:\n\n");
fprintf(fid, "| Metric | Selected | `20260210_180826_modified` | Other / selected |\n");
fprintf(fid, "|---|---:|---:|---:|\n");
write_pair_row(fid, "x1 IAE total", state_iae_schedule(selectedSch, 1), state_iae_schedule(thirdSch, 1));
write_pair_row(fid, "x1 IAE case 1", double(selectedSch.IAE_x1_c1), double(thirdSch.IAE_x1_c1));
write_pair_row(fid, "x1 IAE case 2", double(selectedSch.IAE_x1_c2), double(thirdSch.IAE_x1_c2));
write_pair_row(fid, "x3 IAE case 2", double(selectedSch.IAE_x3_c2), double(thirdSch.IAE_x3_c2));
fprintf(fid, "\nThose ratios are also not exactly 4.4897x. The value `4.4897` is very close to the final-fidelity `IAE_case_c2` value for `20260210_180826` (`%.12g`), but that is a case-total IAE value, not a volume-error ratio.\n\n", ...
    double(thirdF1.noisy_IAE_case_c2));

fprintf(fid, "## Conclusion\n\n");
fprintf(fid, "Selected controller: `%s`.\n\n", selectedScheduleId);
fprintf(fid, "- Better than benchmark in the final f=1 Pareto metrics: `J_track` is %.4f%% lower and `J_TV` is %.4f%% lower.\n", ...
    100 * (1 - selectedFinal.Jtrack / benchmarkFinal.Jtrack), ...
    100 * (1 - selectedFinal.JTV / benchmarkFinal.JTV));
fprintf(fid, "- Best `J_track` position in the main Pareto analysis: it is the lowest-`J_track`/SSE point on the combined run1+run2 non-DOE Pareto frontier.\n");
fprintf(fid, "- Better than benchmark in aggregate schedule TV and runtime: `SSdU_total` is %.4f%% lower and wall time is %.2fx faster.\n", ...
    100 * (1 - double(selectedSch.SSdU_total) / double(benchmarkSch.SSdU_total)), ...
    double(benchmarkSch.wall_time_s) / double(selectedSch.wall_time_s));
fprintf(fid, "- Comparable to benchmark in aggregate schedule tracking: `SSE_total` is %.4f%% higher, which is within 2%%.\n", ...
    100 * (double(selectedSch.SSE_total) / double(benchmarkSch.SSE_total) - 1));
fprintf(fid, "- Better than benchmark in the constant-`Ssp` schedule case for `SSE`, `SSdU`, and biomass tracking: `SSE` is %.4f%% lower, `SSdU` is %.4f%% lower, and `x2` IAE is %.4f%% lower.\n", ...
    100 * (1 - double(selectedSch.SSE_c1) / double(benchmarkSch.SSE_c1)), ...
    100 * (1 - double(selectedSch.SSdU_c1) / double(benchmarkSch.SSdU_c1)), ...
    100 * (1 - double(selectedSch.IAE_x2_c1) / double(benchmarkSch.IAE_x2_c1)));
fprintf(fid, "- Better than benchmark in the moving-`Ssp` schedule case for TV and volume tracking: `SSdU` is %.4f%% lower and `x1` IAE is %.4f%% lower.\n", ...
    100 * (1 - double(selectedSch.SSdU_c2) / double(benchmarkSch.SSdU_c2)), ...
    100 * (1 - double(selectedSch.IAE_x1_c2) / double(benchmarkSch.IAE_x1_c2)));
fprintf(fid, "- Worse than benchmark in some secondary schedule tracking metrics: aggregate `x1` IAE is %.4f%% higher and aggregate `x3` IAE is %.4f%% higher. In the moving-`Ssp` case, `SSE` is %.4f%% higher and `x3` IAE is %.4f%% higher.\n", ...
    100 * (state_iae_schedule(selectedSch, 1) / state_iae_schedule(benchmarkSch, 1) - 1), ...
    100 * (state_iae_schedule(selectedSch, 3) / state_iae_schedule(benchmarkSch, 3) - 1), ...
    100 * (double(selectedSch.SSE_c2) / double(benchmarkSch.SSE_c2) - 1), ...
    100 * (double(selectedSch.IAE_x3_c2) / double(benchmarkSch.IAE_x3_c2) - 1));
fprintf(fid, "- Biomass tracking is the key practical tracking result: aggregate `x2` IAE is only %.4f%% higher than benchmark across the two schedule cases.\n\n", ...
    100 * (state_iae_schedule(selectedSch, 2) / state_iae_schedule(benchmarkSch, 2) - 1));

fprintf(fid, "Alternative controller: `ts_20260210_151703_modified`.\n\n");
fprintf(fid, "- Lower schedule TV than the selected controller: `SSdU_total` is %.4f%% lower.\n", ...
    100 * (1 - double(nextSch.SSdU_total) / double(selectedSch.SSdU_total)));
fprintf(fid, "- Faster than the selected controller: schedule wall time is %.2fx lower/faster.\n", ...
    double(selectedSch.wall_time_s) / double(nextSch.wall_time_s));
fprintf(fid, "- Slightly better aggregate biomass IAE than the selected controller: `x2` IAE is %.4f%% lower.\n", ...
    100 * (1 - state_iae_schedule(nextSch, 2) / state_iae_schedule(selectedSch, 2)));
fprintf(fid, "- Worse aggregate schedule tracking than the selected controller: `SSE_total` is %.4f%% higher.\n", ...
    100 * (double(nextSch.SSE_total) / double(selectedSch.SSE_total) - 1));
fprintf(fid, "- Much worse volume tracking than the selected controller: aggregate `x1` IAE is %.4fx higher; constant-`Ssp` `x1` IAE is %.4fx higher; moving-`Ssp` `x1` IAE is %.4fx higher.\n", ...
    state_iae_schedule(nextSch, 1) / state_iae_schedule(selectedSch, 1), ...
    double(nextSch.IAE_x1_c1) / double(selectedSch.IAE_x1_c1), ...
    double(nextSch.IAE_x1_c2) / double(selectedSch.IAE_x1_c2));
fprintf(fid, "- Much worse substrate tracking than the selected controller: aggregate `x3` IAE is %.4f%% higher; moving-`Ssp` `x3` IAE is %.4f%% higher.\n\n", ...
    100 * (state_iae_schedule(nextSch, 3) / state_iae_schedule(selectedSch, 3) - 1), ...
    100 * (double(nextSch.IAE_x3_c2) / double(selectedSch.IAE_x3_c2) - 1));

fprintf(fid, "Decision rationale.\n\n");
fprintf(fid, "- Prefer the selected controller because the important requirements are low `J_track`, acceptable biomass tracking, reduced TV, and reduced compute time. The selected controller is the lowest-`J_track` Pareto point, is within 2%% schedule `SSE_total` of the benchmark, keeps biomass schedule IAE within %.4f%% of the benchmark, lowers schedule TV by %.4f%% versus benchmark, and is %.2fx faster than benchmark.\n", ...
    100 * (state_iae_schedule(selectedSch, 2) / state_iae_schedule(benchmarkSch, 2) - 1), ...
    100 * (1 - double(selectedSch.SSdU_total) / double(benchmarkSch.SSdU_total)), ...
    double(benchmarkSch.wall_time_s) / double(selectedSch.wall_time_s));
fprintf(fid, "- Do not prefer the adjacent lower-TV alternative because its lower TV and faster runtime come with substantially worse schedule volume and substrate tracking: `x1` IAE is %.4fx higher and `x3` IAE is %.4f%% higher than the selected controller.\n", ...
    state_iae_schedule(nextSch, 1) / state_iae_schedule(selectedSch, 1), ...
    100 * (state_iae_schedule(nextSch, 3) / state_iae_schedule(selectedSch, 3) - 1));
fprintf(fid, "- The volume penalty is less important than biomass tracking here because volume disturbances are not frequent exogenous events in this setup; the in- and outflows are controlled by the NMPC. The substrate penalty is acceptable for the selected controller because biomass remains accurately controlled, and biomass is the primary control objective.\n");
fprintf(fid, "- Therefore, the selected controller is the preferable compromise: it preserves the lowest tracking objective and practical biomass performance while delivering much lower TV and much faster wall time than the benchmark.\n\n");
fprintf(fid, "The old comment block mixed contexts. Current data do not support the exact statements `6.0208%% wall time`, `16.6 times faster`, `8.8700%% lower SSdU`, or `4.4897 times smaller volume tracking error` for the selected controller versus the schedule benchmark or the adjacent lower-`J_TV` Pareto controller.\n");

fprintf("Wrote %s\n", outMd);

function row = require_one(T, mask, label)
idx = find(mask);
if isempty(idx)
    error("Could not find %s.", label);
end
row = T(idx(1), :);
end

function [TallOpt, Tfront] = build_combined_main_front(resultsRoot)
TallOpt = table();
for runKey = ["run1", "run2"]
    T = readtable(fullfile(resultsRoot, runKey, "results.csv"), "TextType", "string");
    T.run_key = repmat(runKey, height(T), 1);
    T.iteration = (1:height(T)).';
    TallOpt = [TallOpt; T(double(T.iteration) > 20, :)]; %#ok<AGROW>
end
isFront = compute_pareto_mask_local(double(TallOpt.SSE), double(TallOpt.SSdU));
TallOpt.is_combined_pareto = isFront;
Tfront = sortrows(TallOpt(isFront, :), "SSE", "ascend");
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

function M = compute_out_metrics(matPath, settlingTol)
S = load(matPath, "out");
out = S.out;
M.Jtrack = double(out.SSE);
M.JTV = double(out.SSdU);
M.runtime_s = nan;
if isfield(out, "runtime_s")
    M.runtime_s = double(out.runtime_s);
end
M.IAE = nan(2, 3);
M.settle_h = nan(2, 3);
for c = 1:min(numel(out.case), 2)
    [settle_h, iae] = summarize_case_metrics_simple(out.case(c), settlingTol);
    M.settle_h(c, :) = settle_h(1:3);
    M.IAE(c, :) = iae(1:3);
end
end

function runtimes = schedule_case_runtimes(matPath)
runtimes = nan(1, 2);
S = load(matPath, "out");
if ~isfield(S, "out") || ~isfield(S.out, "case")
    return
end
for c = 1:min(numel(S.out.case), 2)
    if isfield(S.out.case(c), "runtime_s")
        runtimes(c) = double(S.out.case(c).runtime_s);
    end
end
end

function [settlingTimes_h, IAEByState] = summarize_case_metrics_simple(caseStruct, settlingTol)
Y = double(caseStruct.Y);
Ysp = double(caseStruct.Ysp);
nState = min([size(Y, 2), size(Ysp, 2), 3]);
Y = Y(:, 1:nState);
Ysp = Ysp(:, 1:nState);
if isfield(caseStruct, "dt")
    dt = double(caseStruct.dt);
elseif isfield(caseStruct, "tf")
    dt = double(caseStruct.tf) / max(size(Y, 1) - 1, 1);
else
    dt = 1 / 60;
end
t = (0:size(Y, 1)-1).' * dt;
settlingTimes_h = nan(1, 3);
IAEByState = nan(1, 3);
IAEByState(1:nState) = sum(abs(Y - Ysp), 1) * dt;
for i = 1:nState
    relErr = abs(Y(:, i) - Ysp(:, i)) / max(abs(Ysp(end, i)), 1e-9);
    if relErr(end) <= settlingTol
        idx = find_settling_index(relErr, settlingTol);
        if ~isempty(idx)
            settlingTimes_h(i) = t(idx);
        end
    end
end
end

function idx = find_settling_index(relErr, tol)
idx = [];
for k = 1:numel(relErr)
    if all(relErr(k:end) <= tol)
        idx = k;
        return
    end
end
end

function v = state_iae_schedule(T, stateIdx)
v = double(T.(sprintf("IAE_x%d_c1", stateIdx))) + double(T.(sprintf("IAE_x%d_c2", stateIdx)));
end

function v = state_iae_f1(T, prefix, stateIdx)
v = double(T.(sprintf("%s_IAE_x%d_c1", prefix, stateIdx))) + ...
    double(T.(sprintf("%s_IAE_x%d_c2", prefix, stateIdx)));
end

function write_main_row(fid, label, T)
fprintf(fid, "| %s `%s/%s` | %.12g | %.12g | %.12g | %.12g | %.12g |\n", ...
    label, string(T.run_key), string(T.timestamp), double(T.SSE), double(T.SSdU), double(T.J), double(T.theta_1), double(T.runtime_s));
end

function write_tuning_row(fid, label, T)
theta = theta_vector(T);
Q = 10 .^ theta(4:6);
Ru = 10 .^ theta(7:9);
Rdu = 10 .^ theta(10:12);
fprintf(fid, "| %s `%s` | `%s` | %.12g | %.12g | %.12g | `%s` | `%s` | `%s` | `%s` |\n", ...
    label, string(T.controller_id), string(T.source), theta(1), double(T.m), double(T.p), ...
    fmt_vec(Q), fmt_vec(Ru), fmt_vec(Rdu), fmt_vec(theta));
end

function theta = theta_vector(T)
theta = nan(1, 12);
for k = 1:12
    theta(k) = double(T.(sprintf("theta_%d", k)));
end
end

function s = fmt_vec(v)
parts = strings(1, numel(v));
for k = 1:numel(v)
    if abs(v(k)) < realmin
        parts(k) = "0";
    else
        parts(k) = sprintf("%.6g", v(k));
    end
end
s = "[" + strjoin(parts, ", ") + "]";
end

function write_metric_row(fid, label, selectedValue, benchmarkValue, note)
fprintf(fid, "| %s | %.12g | %.12g | %.6g | %s (%s) |\n", ...
    label, selectedValue, benchmarkValue, selectedValue / benchmarkValue, relation_text(selectedValue, benchmarkValue), note);
end

function write_case_metric_row(fid, caseIdx, caseLabel, label, selectedValue, benchmarkValue)
fprintf(fid, "| %s | %s | %s | %.12g | %.12g | %.6g |\n", string(caseIdx), caseLabel, label, selectedValue, benchmarkValue, selectedValue / benchmarkValue);
end

function label = schedule_case_label(caseIdx)
if caseIdx == 1
    label = "constant/steady-state Ssp, upper panels";
elseif caseIdx == 2
    label = "moving Ssp, lower panels";
else
    label = "unknown";
end
end

function write_pair_row(fid, label, selectedValue, otherValue)
fprintf(fid, "| %s | %.12g | %.12g | %.6g |\n", label, selectedValue, otherValue, otherValue / selectedValue);
end

function txt = relation_text(value, reference)
ratio = value / reference;
if ratio <= 1
    txt = sprintf("%s lower", pct_delta(ratio));
else
    txt = sprintf("%s higher", pct_delta(ratio));
end
end

function txt = pct_higher(value, reference)
txt = sprintf("%.4f%% higher", 100 * (value / reference - 1));
end

function txt = pct_lower(value, reference)
txt = sprintf("%.4f%% lower", 100 * (1 - value / reference));
end

function txt = pct_of(value, reference)
txt = sprintf("%.4f%%", 100 * value / reference);
end

function txt = pct_delta(ratio)
txt = sprintf("%.4f%%", 100 * abs(ratio - 1));
end
