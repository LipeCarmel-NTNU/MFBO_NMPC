%% verify_paper_values
% Recomputes every derived number reported in Overleaf---MFBO-NMPC/main.tex
% directly from the stored result files, and prints them as ground truth.
%
% Companion document: Overleaf---MFBO-NMPC/value_dictionary.md, which lists
% each value in the manuscript together with the source it comes from. Every
% quantity printed here is tagged with the manuscript location that quotes it.
%
% Run from anywhere; paths are resolved relative to this file.

clear; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
res  = fullfile(root, 'results');

fprintf('Repository root: %s\n\n', root);

%% Load the two BO runs
% results.csv columns: timestamp, SSE, SSdU, J, runtime_s, theta_1..theta_12
% theta = [z, theta_p, theta_m, q1..q3, r_u1..r_u3, r_du1..r_du3]
N_INIT = 20;                      % DOE evaluations per run (main.py: N_INIT)

cases = struct('label', {'Case 1', 'Case 2'}, 'dir', {'run1', 'run2'});
R = struct([]);
for c = 1:numel(cases)
    T = readtable(fullfile(res, cases(c).dir, 'results.csv'), ...
                  'TextType', 'string', 'VariableNamingRule', 'preserve');
    n = height(T);
    R(c).label   = cases(c).label;
    R(c).ts      = string(T.timestamp);
    R(c).SSE     = T.SSE;
    R(c).SSdU    = T.SSdU;
    R(c).t_s     = T.runtime_s;
    R(c).z       = T.theta_1;
    R(c).Nc      = round(T.theta_3) + 1;
    R(c).Np      = round(T.theta_2) + R(c).Nc;
    R(c).isDOE   = (1:n)' <= N_INIT;
    R(c).theta   = table2array(T(:, 6:17));
    fprintf('Loaded %s: %d evaluations (%d DOE + %d optimisation)\n', ...
            R(c).label, n, sum(R(c).isDOE), sum(~R(c).isDOE));
end
fprintf('\n');

%% Runtime and fidelity allocation
% Quoted in: Results 4.1 ("Computational cost and fidelity allocation"),
% Discussion "How runtime awareness buys efficiency", discussion.tex
% "Computational allocation and fidelity deployment".
fprintf('== Runtime and fidelity allocation ==\n');
doe_all = 0; opt_all = 0; z_opt_all = [];
for c = 1:2
    d = sum(R(c).t_s( R(c).isDOE)) / 60;
    o = sum(R(c).t_s(~R(c).isDOE)) / 60;
    T = d + o;
    doe_all = doe_all + d;  opt_all = opt_all + o;
    z_opt_all = [z_opt_all; R(c).z(~R(c).isDOE)];               %#ok<AGROW>

    [tmax, imax] = max(R(c).t_s);
    tmax_opt = max(R(c).t_s(~R(c).isDOE));

    fprintf('%s\n', R(c).label);
    fprintf('  DOE            = %9.3f min (%6.3f %%)\n', d, 100*d/T);
    fprintf('  optimisation   = %9.3f min (%6.3f %%)\n', o, 100*o/T);
    fprintf('  total          = %9.3f min = %6.4f h\n', T, T/60);
    fprintf('  mean per eval  : DOE %7.3f min | optimisation %7.3f min\n', ...
            d/sum(R(c).isDOE), o/sum(~R(c).isDOE));
    fprintf('  most expensive evaluation : %6.4f h = %6.3f %% of run total (%s, idx %d)\n', ...
            tmax/3600, 100*(tmax/60)/T, ternary(R(c).isDOE(imax), 'DOE', 'optimisation'), imax);
    fprintf('  slowest optimisation-phase evaluation : %7.3f min\n', tmax_opt/60);
    fprintf('  mean z : DOE %8.6f | optimisation %8.6f | all %8.6f\n', ...
            mean(R(c).z(R(c).isDOE)), mean(R(c).z(~R(c).isDOE)), mean(R(c).z));
    fprintf('  z range over all evaluations : [%8.6f, %8.6f]\n', min(R(c).z), max(R(c).z));
    nc1 = sum(R(c).Nc(~R(c).isDOE) == 1);
    fprintf('  N_c = 1 among optimisation evaluations : %d/%d (%7.4f %%)\n\n', ...
            nc1, sum(~R(c).isDOE), 100*nc1/sum(~R(c).isDOE));
end
tot_all = doe_all + opt_all;
zbar = mean(z_opt_all);
fprintf('Combined\n');
fprintf('  DOE          = %9.3f min (%6.4f %%)\n', doe_all, 100*doe_all/tot_all);
fprintf('  optimisation = %9.3f min (%6.4f %%)\n', opt_all, 100*opt_all/tot_all);
fprintf('  total        = %9.3f min\n', tot_all);
fprintf('  mean per eval: DOE %7.3f min | optimisation %7.3f min\n', doe_all/40, opt_all/202);
fprintf('  mean optimisation-phase z (zbar) = %8.6f  ->  zbar*T = %6.4f h (T = 10 h)\n', zbar, 10*zbar);
r_cost = (sum(R(2).t_s(~R(2).isDOE))/101) / (sum(R(1).t_s(~R(1).isDOE))/101);
fprintf('  Case 2 / Case 1 mean optimisation-phase cost ratio = %6.4f\n\n', r_cost);

%% Pareto sets
% Four distinct non-dominated sets circulate in the manuscript. They are all
% computed here so that each reported count can be attributed unambiguously.
fprintf('== Pareto sets ==\n');
allSSE  = [R(1).SSE;  R(2).SSE];
allSSdU = [R(1).SSdU; R(2).SSdU];
allNc   = [R(1).Nc;   R(2).Nc];
allNp   = [R(1).Np;   R(2).Np];
allz    = [R(1).z;    R(2).z];
allTs   = [R(1).ts;   R(2).ts];
allCase = [repmat("Case 1", numel(R(1).ts), 1); repmat("Case 2", numel(R(2).ts), 1)];
allDOE  = [R(1).isDOE; R(2).isDOE];

% (i) merged over all 242 evaluations, DOE included: basis of the N_c/N_p boxplots
sel = pareto_min([allSSE, allSSdU]);
report_pareto('merged, all 242 evaluations (N_c/N_p boxplots, Results 4.2)', ...
              allCase(sel), allTs(sel), allz(sel), allNc(sel), allNp(sel), ...
              allSSE(sel), allSSdU(sel));

% (ii) merged over the 202 optimisation-phase evaluations: the set the
% Algorithm-summary subsection defines, and the validation-schedule cohort
keep = ~allDOE;
sub  = find(keep);
sel2 = sub(pareto_min([allSSE(keep), allSSdU(keep)]));
report_pareto('merged, optimisation phase only (final Pareto set, validation cohort)', ...
              allCase(sel2), allTs(sel2), allz(sel2), allNc(sel2), allNp(sel2), ...
              allSSE(sel2), allSSdU(sel2));

% (iii)/(iv) within-run sets
for c = 1:2
    s = pareto_min([R(c).SSE, R(c).SSdU]);
    fprintf('  within-run %s, all 121 evaluations : %d members (N_c = 1 in %d)\n', ...
            R(c).label, sum(s), sum(R(c).Nc(s) == 1));
    o = ~R(c).isDOE;  io = find(o);
    s2 = io(pareto_min([R(c).SSE(o), R(c).SSdU(o)]));
    fprintf('  within-run %s, optimisation phase only : %d members (N_c = 1 in %d, N_p = 1 in %d)\n', ...
            R(c).label, numel(s2), sum(R(c).Nc(s2) == 1), sum(R(c).Np(s2) == 1));
end
fprintf('\n');

%% Benchmark comparison at full fidelity
% Quoted in: Results 4.3, Table 1, Abstract, Introduction, Conclusions.
fprintf('== Benchmark comparison (full fidelity, same noise realisation) ==\n');
B = readtable(fullfile(res, 'benchmark_reference_controller', ...
              'benchmark_full_f1_same_noise_fix', 'results_benchmark.csv'), ...
              'TextType', 'string', 'VariableNamingRule', 'preserve');
bSSE  = B.SSE;
bSSdU = B.SSdU;
bNc   = round(B.theta_3) + 1;
bNp   = round(B.theta_2) + bNc;
fprintf('Benchmark: J_track = %10.4f | J_du = %10.6f | N_c = %d | N_p = %d\n', ...
        bSSE, bSSdU, bNc, bNp);
fprintf('  Q   = [%g %g %g]   (rescaled to Q11 = 1: [%g %g %g])\n', ...
        10.^B.theta_4, 10.^B.theta_5, 10.^B.theta_6, ...
        1, 10.^(B.theta_5-B.theta_4), 10.^(B.theta_6-B.theta_4));
fprintf('  Ru  = [%g %g %g]\n', 10.^B.theta_7, 10.^B.theta_8, 10.^B.theta_9);
fprintf('  Rdu = [%g %g %g]   (rescaled: [%g %g %g])\n\n', ...
        10.^B.theta_10, 10.^B.theta_11, 10.^B.theta_12, ...
        10.^(B.theta_10-B.theta_4), 10.^(B.theta_11-B.theta_4), 10.^(B.theta_12-B.theta_4));

for c = 1:2
    F = readtable(fullfile(res, 'final_fidelity_same_noise', ...
                  sprintf('run%d_full_f1_same_noise', c), 'results_full.csv'), ...
                  'TextType', 'string', 'VariableNamingRule', 'preserve');
    F = sortrows(F, 'SSE');
    fprintf('%s, refined at z = 1:\n', R(c).label);
    for i = 1:height(F)
        nc = round(F.theta_3(i)) + 1;
        np = round(F.theta_2(i)) + nc;
        fprintf(['  %s  J_track = %10.3f  J_du = %9.6f  ratio = %8.3f  ' ...
                 'N_c = %2d  N_p = %2d  Q = [%.4g %.4g %.4g]  Rdu = [%.4g %.4g %.4g]\n'], ...
                F.timestamp(i), F.SSE(i), F.SSdU(i), bSSdU/F.SSdU(i), nc, np, ...
                10^F.theta_4(i), 10^F.theta_5(i), 10^F.theta_6(i), ...
                10^F.theta_10(i), 10^F.theta_11(i), 10^F.theta_12(i));
    end
    if c == 2
        m = F.SSE < 12300;                     % the two "essentially identical tracking" members
        fprintf('  headline range: J_track in [%.3f, %.3f], J_du in [%.6f, %.6f]\n', ...
                min(F.SSE(m)), max(F.SSE(m)), min(F.SSdU(m)), max(F.SSdU(m)));
        fprintf('  control-variation ratio vs benchmark: %.4f to %.4f  ->  reported as "%d to %d times"\n', ...
                bSSdU/max(F.SSdU(m)), bSSdU/min(F.SSdU(m)), ...
                round(bSSdU/max(F.SSdU(m))), round(bSSdU/min(F.SSdU(m))));
    end
    fprintf('\n');
end

%% Settling analysis
% Quoted in: Results 4.4, Figure "python_settling_time_boxplot".
% Cohort: 21 controllers x 2 initial-condition scenarios = 42 evaluations,
% re-evaluated at z = 1 without measurement noise.
fprintf('== Settling analysis (reference horizon zbar*T) ==\n');
S = readtable(fullfile(here, 'boxplot_settling_time_data.csv'), ...
              'TextType', 'string', 'VariableNamingRule', 'preserve');
ref_h = 10 * zbar;
key = S.run_label + "|" + S.timestamp + "|" + string(S.case_id);
units = unique(key, 'stable');
states = ["x1", "x2", "x3"];
settled = false(numel(units), 3);
isnanv  = false(numel(units), 3);
val = S.value;
if ~isnumeric(val); val = str2double(string(val)); end     % 'NaN' entries may import as text
for u = 1:numel(units)
    for s = 1:3
        v = val(key == units(u) & S.state == states(s));
        isnanv(u,s)  = isempty(v) || isnan(v(1));
        settled(u,s) = ~isnanv(u,s) && v(1) <= ref_h;
    end
end
n = numel(units);
fprintf('Reference horizon = %.5f h; cohort = %d controller-scenario pairs\n', ref_h, n);
lbl = ["V", "X", "S"];
for s = 1:3
    fprintf('  %s: settled before reference : %5.2f %% (%d/%d) | never within 5%% band : %5.2f %% (%d/%d)\n', ...
            lbl(s), 100*sum(settled(:,s))/n, sum(settled(:,s)), n, ...
            100*sum(isnanv(:,s))/n, sum(isnanv(:,s)), n);
end
fprintf('  all three states simultaneously : %5.2f %% (%d/%d)\n\n', ...
        100*sum(all(settled,2))/n, sum(all(settled,2)), n);

%% Extrapolation validation
% Quoted in: Methods "Multi-fidelity closed-loop evaluation" (R^2) and
% Results 4.4 (endpoint errors, coverage).
fprintf('== Extrapolation validation ==\n');
coef_file = fullfile(res, 'numerical results', 'surrogate_cheb_coeffs.txt');
disp(strjoin(readlines(coef_file), newline));

L = readlines(fullfile(res, 'numerical results', 'surrogate_test_summary.txt'));
h = find(startsWith(L, 'run,file,'), 1);
rows = L(h+1:end);
rows = rows(strlength(strip(rows)) > 0);
parts = split(rows, ',');
z_eval  = str2double(parts(:,3));
err_dU  = abs(1 - str2double(parts(:,4)));
err_E   = abs(1 - str2double(parts(:,5)));
fprintf('Processed files: %d (z in [%.6f, %.6f], mean %.6f)\n', ...
        numel(z_eval), min(z_eval), max(z_eval), mean(z_eval));
fprintf('  J_du   endpoint |error|: mean %10.4e  median %10.4e  p95 %10.4e  max %10.4e\n', ...
        mean(err_dU), median(err_dU), pctl(err_dU, 0.95), max(err_dU));
fprintf('  J_track endpoint |error|: mean %10.4e  median %10.4e  p95 %10.4e  max %10.4e\n', ...
        mean(err_E), median(err_E), pctl(err_E, 0.95), max(err_E));
fprintf('  within 1 %% for both objectives: %5.2f %% (%d/%d)\n\n', ...
        100*mean(err_dU <= 0.01 & err_E <= 0.01), ...
        sum(err_dU <= 0.01 & err_E <= 0.01), numel(z_eval));

%% Validation schedule
% Quoted in: Results 4.5 and the parked Draft Helper item on volume IAE.
fprintf('== Validation schedule (30 h, Xsp 7 -> 13 -> 16) ==\n');
V = readtable(fullfile(res, 'numerical results', 'setpoint_schedule_metrics_same_noise.csv'), ...
              'TextType', 'string', 'VariableNamingRule', 'preserve');
isB = V.is_benchmark == 1;
% Policy (a) = steady-state substrate policy -> columns *_c1
% Policy (b) = Ssp(X) feedback policy        -> columns *_c2
for pol = 1:2
    col = sprintf('SSdU_c%d', pol);
    v = V.(col);
    [vs, ord] = sort(v, 'descend');
    fprintf('Policy (%c): highest = %s (%.6f), second = %s (%.6f), ratio = %.4f\n', ...
            char('a' + pol - 1), V.controller_id(ord(1)), vs(1), ...
            V.controller_id(ord(2)), vs(2), vs(1)/vs(2));
end
IAE_V = V.IAE_x1_c1 + V.IAE_x1_c2;
[iv, ord] = sort(IAE_V);
fprintf('Volume IAE summed over both policies (Pareto controllers only):\n');
for i = 1:min(3, numel(ord))
    if isB(ord(i)); continue; end
    fprintf('  %-34s %.4f\n', V.controller_id(ord(i)), iv(i));
end
nonB = find(~isB);
ivn = sort(IAE_V(nonB));
fprintf('  best = %.4f, second best = %.4f, factor = %.3f\n', ivn(1), ivn(2), ivn(2)/ivn(1));
fprintf('  benchmark for reference = %.4f\n', IAE_V(isB));

% Biomass settling lag after the final setpoint step at t = 20 h
lag_a = V.settle_x2_h_c1 - 20;
lag_b = V.settle_x2_h_c2 - 20;
for pol = 1:2
    lag = ternary(pol == 1, lag_a, lag_b);
    bad = ~isB & (isnan(lag) | lag > 3);
    fprintf('Policy (%c): %d of %d Pareto controllers took > 3 h to settle X or never settled\n', ...
            char('a' + pol - 1), sum(bad), sum(~isB));
end
fprintf('\nDone.\n');

%% Local functions
function idx = pareto_min(F)
    % Logical index of the non-dominated rows of F under simultaneous minimisation.
    n = size(F, 1);
    idx = true(n, 1);
    for i = 1:n
        others = [1:i-1, i+1:n];
        dominated = all(F(others,:) <= F(i,:), 2) & any(F(others,:) < F(i,:), 2);
        idx(i) = ~any(dominated);
    end
end

function report_pareto(name, cs, ts, z, Nc, Np, SSE, SSdU)
    [~, ord] = sort(SSdU);
    fprintf('  %s: %d members (Case 1: %d, Case 2: %d)\n', ...
            name, numel(ord), sum(cs == "Case 1"), sum(cs == "Case 2"));
    for i = ord'
        fprintf('    %s %s  z = %.4f  N_c = %2d  N_p = %2d  J_track = %10.3f  J_du = %.6f\n', ...
                cs(i), ts(i), z(i), Nc(i), Np(i), SSE(i), SSdU(i));
    end
    for c = ["Case 1", "Case 2"]
        m = cs == c;
        if ~any(m); continue; end
        fprintf('    %s: N_c = %s, N_p = %s | medians N_c = %g, N_p = %g | N_c = 1 in %d/%d\n', ...
                c, mat2str(Nc(m)'), mat2str(Np(m)'), median(Nc(m)), median(Np(m)), ...
                sum(Nc(m) == 1), sum(m));
    end
    fprintf('    combined N_c = 1 in %d/%d | min z = %.4f\n', ...
            sum(Nc == 1), numel(Nc), min(z));
    fprintf('    Case 1 min J_du = %.6f, max J_track = %.3f | Case 2 min J_track = %.3f, J_du in [%.4f, %.4f]\n\n', ...
            min(SSdU(cs == "Case 1")), max(SSE(cs == "Case 1")), ...
            min(SSE(cs == "Case 2")), min(SSdU(cs == "Case 2")), max(SSdU(cs == "Case 2")));
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

function q = pctl(v, p)
    % Linear-interpolation percentile, matching the convention used to produce
    % the p95 figures quoted in the manuscript. Avoids a toolbox dependency.
    v = sort(v(:));
    k = (numel(v) - 1) * p;
    lo = floor(k) + 1;
    hi = min(lo + 1, numel(v));
    q = v(lo) + (v(hi) - v(lo)) * (k - floor(k));
end
