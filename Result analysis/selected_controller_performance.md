# Selected controller performance audit

Selected controller: `ts_20260211_122653_modified` in the schedule validation tables. Its base timestamp is `20260211_122653`, from `run2`.

The same timestamp appears in `main_pareto` as `run2/20260211_122653`. These are different evaluation contexts and should not be mixed:

- Pareto optimization / final-fidelity validation: original run1/run2 objective data, then f=1 final-fidelity re-evaluation against `benchmark_full_f1_same_noise_fix`.
- Schedule validation: setpoint schedule `setpoint_schedule_xsp_7_13_16/same_noise`, using the `_modified` controller row and `benchmark_fix` from the schedule table. The schedule analysis plots two schedule cases: case 1 is the upper-row constant/steady-state `Ssp` case, and case 2 is the lower-row moving `Ssp` case.

The selected controller is a good choice because it dominates the benchmark in the final f=1 Pareto metrics: `J_track` is 0.3139% lower and `J_TV` is 88.9916% lower than the benchmark. In schedule validation, the same controller is 14.25x faster than the schedule benchmark, `SSE_total` is within 2% of the benchmark, `SSdU_total` is 31.6723% lower, and biomass (`x2`) IAE is only 0.2246% higher. The adjacent lower-`J_TV` schedule controller is tempting because it gives lower `SSdU_total`, but its schedule volume (`x1`) IAE is 2.8893x higher and substrate (`x3`) IAE is 51.4328% higher than the selected controller.

## Tuning parameters

The schedule/final-fidelity validation uses the selected controller at `z = 1`. The original BO point in `main_pareto` had `z = 0.842711112439` before final-fidelity re-evaluation.

| Controller | source | z/f | m | p | Q diag | R_u diag | R_du diag | theta |
|---|---|---:|---:|---:|---|---|---|---|
| selected `ts_20260211_122653_modified` | `run2` | 1 | 1 | 15 | `[1, 0.00229222, 0.0453124]` | `[0, 0, 0]` | `[0.0193382, 0.14944, 0.0022462]` | `[1, 14, 0, 0, -2.63974, -1.34378, -1000, -1000, -1000, -1.71358, -0.825533, -2.64855]` |
| benchmark `benchmark_fix` | `benchmark_fix` | 1 | 6 | 61 | `[10, 1, 1]` | `[0, 0, 0]` | `[10, 10, 10]` | `[1, 55, 5, 1, 0, 0, -1000, -1000, -1000, 1, 1, 1]` |

## Pareto optimization context

In the combined non-DOE `main_pareto` frontier, `20260211_122653` is the lowest `J_track` point and the highest `J_TV` point on the frontier.

| Scope | Rank by lowest J_track/SSE | Rank by lowest J_TV/SSdU |
|---|---:|---:|
| run2 optimization rows | 1 of 101 | 35 of 101 |
| run2 Pareto frontier | 1 of 10 | 10 of 10 |
| combined run1+run2 optimization rows | 1 of 202 | 76 of 202 |
| combined main Pareto frontier | 1 of 10 | 10 of 10 |

Original optimization objective values, before f=1 final-fidelity re-evaluation:

| Controller | J_track/SSE | J_TV/SSdU | J | z/f | Runtime s |
|---|---:|---:|---:|---:|---:|
| selected `run2/20260211_122653` | 12214.913069 | 0.368936444696 | 15904.277516 | 0.842711112439 | 357.13901 |
| next lower J_TV frontier point `run2/20260210_151703` | 12273.2069237 | 0.178714926414 | 14060.3561879 | 0.801429873974 | 331.034231 |

Moving from the selected controller to the next lower-`J_TV` Pareto point (`20260210_151703`) would reduce original `J_TV` by 51.5594%, but increase original `J_track` by 0.4772%.

Final f=1 noisy validation against the benchmark:

| Metric | Selected | Benchmark | Selected / benchmark | Interpretation |
|---|---:|---:|---:|---|
| J_track/SSE | 12220.9890421 | 12259.4702149 | 0.996861 | 0.3139% lower (lower is better) |
| J_TV/SSdU | 0.387029255072 | 3.51575951978 | 0.110084 | 88.9916% lower (lower is better) |
| runtime_s | 434.116606 | 8892.896247 | 0.0488161 | 95.1184% lower (lower is better) |
| x1 IAE total | 0.203381042567 | 0.191316107644 | 1.06306 | 6.3063% higher (lower is better) |
| x2 IAE total | 31.8306365661 | 31.3422805535 | 1.01558 | 1.5581% higher (lower is better) |
| x3 IAE total | 1.75976257035 | 2.23048280735 | 0.78896 | 21.1040% lower (lower is better) |

The final f=1 comparison gives `J_track = 12220.9890421`, which is 0.3139% lower than the benchmark, and `J_TV = 0.387029255072`, which is 88.9916% lower than the benchmark.

Adjacent final f=1 Pareto point check (`20260210_151703`):

| Metric | Selected `20260211_122653` | Next `20260210_151703` | Next / selected |
|---|---:|---:|---:|
| noisy J_track | 12220.9890421 | 12298.7816848 | 1.00637 |
| noisy J_TV | 0.387029255072 | 0.162122479332 | 0.418889 |
| x1 IAE total | 0.203381042567 | 0.428370012771 | 2.10624 |
| x2 IAE total | 31.8306365661 | 32.7809699261 | 1.02986 |
| x3 IAE total | 1.75976257035 | 2.44883180364 | 1.39157 |

For the adjacent final-fidelity Pareto point, the selected controller has 2.10624x smaller volume IAE than `20260210_151703`, not 4.4897x.

## Schedule validation context

Schedule validation uses `ts_20260211_122653_modified`, from `C:\Users\lipec\OneDrive - NTNU\Documents\GitHub\MFBO_NMPC\results\numerical results\setpoint_schedule_metrics_same_noise.csv`, and benchmark row `benchmark_fix`.

The two rows in `analyze_setpoint_schedule_metrics.m` are two schedule cases from each `out.case` array, not two separate result folders. Case 1, plotted in the upper panels, keeps the steady-state substrate setpoint returned by the setpoint/terminal-cost workflow. Case 2, plotted in the lower panels, overwrites substrate setpoint online using `Ssp = min(3, 2*(Xsp - X))`, so `Ssp` moves during the trajectory.

The headline schedule numbers in the introduction and conclusion use the summed aggregate over both cases. The average rows below are arithmetic means of case 1 and case 2, included only to show the typical per-case scale.

Aggregate over both schedule cases:

| Metric | Selected | Benchmark | Selected / benchmark | Interpretation |
|---|---:|---:|---:|---|
| SSE_total | 14490.9659118 | 14277.8597015 | 1.01493 | 1.4926% higher (lower is better) |
| SSdU_total | 4.77237861422 | 6.98454927859 | 0.683277 | 31.6723% lower (lower is better) |
| wall_time_s | 2466.02493 | 35132.413512 | 0.0701923 | 92.9808% lower (lower is better) |
| x1 IAE total | 0.162887142157 | 0.150199659224 | 1.08447 | 8.4471% higher (lower is better) |
| x2 IAE total | 63.5176856374 | 63.3753326589 | 1.00225 | 0.2246% higher (lower is better) |
| x3 IAE total | 8.65761665455 | 7.72685926547 | 1.12046 | 12.0457% higher (lower is better) |

Schedule validation case-by-case metrics:

| Case | Ssp policy / plot row | Metric | Selected | Benchmark | Selected / benchmark |
|---:|---|---|---:|---:|---:|
| 1 | constant/steady-state Ssp, upper panels | SSE | 6823.78224428 | 6848.3541805 | 0.996412 |
| 1 | constant/steady-state Ssp, upper panels | SSdU | 1.83502783424 | 3.6926355396 | 0.496943 |
| 1 | constant/steady-state Ssp, upper panels | runtime_s | 803.8258 | 25607.836824 | 0.0313898 |
| 1 | constant/steady-state Ssp, upper panels | x1 IAE | 0.0709553716584 | 0.0366736873808 | 1.93478 |
| 1 | constant/steady-state Ssp, upper panels | x2 IAE | 30.6877445309 | 30.97395953 | 0.990759 |
| 1 | constant/steady-state Ssp, upper panels | x3 IAE | 3.5446210108 | 3.337062955 | 1.0622 |
| 2 | moving Ssp, lower panels | SSE | 7667.18366754 | 7429.50552099 | 1.03199 |
| 2 | moving Ssp, lower panels | SSdU | 2.93735077998 | 3.29191373899 | 0.892293 |
| 2 | moving Ssp, lower panels | runtime_s | 1662.19913 | 9524.576688 | 0.174517 |
| 2 | moving Ssp, lower panels | x1 IAE | 0.0919317704988 | 0.113525971843 | 0.809786 |
| 2 | moving Ssp, lower panels | x2 IAE | 32.8299411065 | 32.4013731289 | 1.01323 |
| 2 | moving Ssp, lower panels | x3 IAE | 5.11299564375 | 4.38979631047 | 1.16475 |
| avg | average of constant and moving Ssp cases | SSE | 7245.48295591 | 7138.92985074 | 1.01493 |
| avg | average of constant and moving Ssp cases | SSdU | 2.38618930711 | 3.49227463929 | 0.683277 |
| avg | average of constant and moving Ssp cases | runtime_s | 1233.012465 | 17566.206756 | 0.0701923 |
| avg | average of constant and moving Ssp cases | x1 IAE | 0.0814435710786 | 0.0750998296119 | 1.08447 |
| avg | average of constant and moving Ssp cases | x2 IAE | 31.7588428187 | 31.6876663295 | 1.00225 |
| avg | average of constant and moving Ssp cases | x3 IAE | 4.32880832728 | 3.86342963274 | 1.12046 |

The schedule validation benchmark comparison is therefore: SSE is 1.4926% higher, SSdU is 31.6723% lower, wall time is 7.0192% of benchmark (14.25x faster), and biomass IAE is 0.2246% higher.

## Check of the `4.4897 times smaller volume tracking error` note

I checked the likely interpretation that this compares the selected controller with the next Pareto controller if we try to lower `J_TV` a little.

Schedule validation, selected vs adjacent lower-`J_TV` schedule controller (`20260210_151703_modified`):

| Metric | Selected | Adjacent lower-J_TV controller | Adjacent / selected |
|---|---:|---:|---:|
| SSE_total | 14490.9659118 | 14785.7261945 | 1.02034 |
| SSdU_total | 4.77237861422 | 1.55993935719 | 0.326868 |
| wall_time_s | 2466.02493 | 1843.84414 | 0.747699 |
| x1 IAE total | 0.162887142157 | 0.470627028825 | 2.88928 |
| x2 IAE total | 63.5176856374 | 63.2473578256 | 0.995744 |
| x3 IAE total | 8.65761665455 | 13.1104740205 | 1.51433 |

Under this interpretation, the selected controller has 2.88928x smaller volume IAE than the adjacent lower-`J_TV` controller. This does not reproduce 4.4897x.

The closest nearby schedule comparison is against `20260210_180826_modified`, not the immediate adjacent lower-`J_TV` point:

| Metric | Selected | `20260210_180826_modified` | Other / selected |
|---|---:|---:|---:|
| x1 IAE total | 0.162887142157 | 0.790851225996 | 4.85521 |
| x1 IAE case 1 | 0.0709553716584 | 0.37081037186 | 5.22597 |
| x1 IAE case 2 | 0.0919317704988 | 0.420040854136 | 4.56905 |
| x3 IAE case 2 | 5.11299564375 | 3.97228954192 | 0.776901 |

Those ratios are also not exactly 4.4897x. The value `4.4897` is very close to the final-fidelity `IAE_case_c2` value for `20260210_180826` (`4.49095063478`), but that is a case-total IAE value, not a volume-error ratio.

## Conclusion

Selected controller: `ts_20260211_122653_modified`.

- Better than benchmark in the final f=1 Pareto metrics: `J_track` is 0.3139% lower and `J_TV` is 88.9916% lower.
- Best `J_track` position in the main Pareto analysis: it is the lowest-`J_track`/SSE point on the combined run1+run2 non-DOE Pareto frontier.
- Better than benchmark in aggregate schedule TV and runtime: `SSdU_total` is 31.6723% lower and wall time is 14.25x faster.
- Comparable to benchmark in aggregate schedule tracking: `SSE_total` is 1.4926% higher, which is within 2%.
- Better than benchmark in the constant-`Ssp` schedule case for `SSE`, `SSdU`, and biomass tracking: `SSE` is 0.3588% lower, `SSdU` is 50.3057% lower, and `x2` IAE is 0.9241% lower.
- Better than benchmark in the moving-`Ssp` schedule case for TV and volume tracking: `SSdU` is 10.7707% lower and `x1` IAE is 19.0214% lower.
- Worse than benchmark in some secondary schedule tracking metrics: aggregate `x1` IAE is 8.4471% higher and aggregate `x3` IAE is 12.0457% higher. In the moving-`Ssp` case, `SSE` is 3.1991% higher and `x3` IAE is 16.4746% higher.
- Biomass tracking is the key practical tracking result: aggregate `x2` IAE is only 0.2246% higher than benchmark across the two schedule cases.

Alternative controller: `ts_20260210_151703_modified`.

- Lower schedule TV than the selected controller: `SSdU_total` is 67.3132% lower.
- Faster than the selected controller: schedule wall time is 1.34x lower/faster.
- Slightly better aggregate biomass IAE than the selected controller: `x2` IAE is 0.4256% lower.
- Worse aggregate schedule tracking than the selected controller: `SSE_total` is 2.0341% higher.
- Much worse volume tracking than the selected controller: aggregate `x1` IAE is 2.8893x higher; constant-`Ssp` `x1` IAE is 3.0384x higher; moving-`Ssp` `x1` IAE is 2.7742x higher.
- Much worse substrate tracking than the selected controller: aggregate `x3` IAE is 51.4328% higher; moving-`Ssp` `x3` IAE is 49.4030% higher.

Decision rationale.

- Prefer the selected controller because the important requirements are low `J_track`, acceptable biomass tracking, reduced TV, and reduced compute time. The selected controller is the lowest-`J_track` Pareto point, is within 2% schedule `SSE_total` of the benchmark, keeps biomass schedule IAE within 0.2246% of the benchmark, lowers schedule TV by 31.6723% versus benchmark, and is 14.25x faster than benchmark.
- Do not prefer the adjacent lower-TV alternative because its lower TV and faster runtime come with substantially worse schedule volume and substrate tracking: `x1` IAE is 2.8893x higher and `x3` IAE is 51.4328% higher than the selected controller.
- The volume penalty is less important than biomass tracking here because volume disturbances are not frequent exogenous events in this setup; the in- and outflows are controlled by the NMPC. The substrate penalty is acceptable for the selected controller because biomass remains accurately controlled, and biomass is the primary control objective.
- Therefore, the selected controller is the preferable compromise: it preserves the lowest tracking objective and practical biomass performance while delivering much lower TV and much faster wall time than the benchmark.

The old comment block mixed contexts. Current data do not support the exact statements `6.0208% wall time`, `16.6 times faster`, `8.8700% lower SSdU`, or `4.4897 times smaller volume tracking error` for the selected controller versus the schedule benchmark or the adjacent lower-`J_TV` Pareto controller.
