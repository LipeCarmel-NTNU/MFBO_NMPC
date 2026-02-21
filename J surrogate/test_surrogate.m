% test_surrogate.m
% Validate Chebyshev-based multifidelity extrapolation against stored run outputs.
%
% What this script checks:
% 1) Time-varying projected full-horizon costs from partial sums:
%       J_hat(t) = J_partial(0:t) / frac(t)
% 2) Normalized prediction trajectories:
%       J_hat(t) / J_final
%    (flat at 1 is ideal)
% 3) End-of-run bias statistics and bias vs fidelity.

clear; close all; clc

scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
fontSize = 16;
plotColors = good_colors(4);

% ------------------------------------------------------------
% Configuration
% ------------------------------------------------------------
fullHorizonHours = 10.0;            % assumed in main setup: T = 10*z hours
maxFilesPerRun = 10;                % keep plots readable (set [] for all)
useCoeffFileIfPresent = false;

    cSSdU = [ ...
         6.442657e-01 ...
         4.682368e-01 ...
        -1.455242e-01 ...
         2.457724e-02 ...
         2.198219e-02 ...
        -1.598959e-02 ...
    ].';

    cSSE = [ ...
         7.827330e-01 ...
         3.771709e-01 ...
        -2.433549e-01 ...
         1.091471e-01 ...
        -2.846259e-02 ...
         1.358572e-03 ...
    ].';

if useCoeffFileIfPresent
    coeffTxt = fullfile(projectRoot, "results", "numerical results", "surrogate_cheb_coeffs.txt");
    [ok, cSSdU_txt, cSSE_txt] = try_parse_coeff_txt(coeffTxt);
    if ok
        cSSdU = cSSdU_txt;
        cSSE = cSSE_txt;
    end
end

% ------------------------------------------------------------
% Load run files
% ------------------------------------------------------------
runDirs = {
    fullfile(projectRoot, "results", "run1"), "Case 1";
    fullfile(projectRoot, "results", "run2"), "Case 2";
};

records = struct( ...
    "runLabel", {}, "fileName", {}, "path", {}, ...
    "tauSSdU", {}, "tauSSE", {}, ...
    "zEval", {}, "ratioSSdU", {}, "ratioSSE", {}, ...
    "endRatioSSdU", {}, "endRatioSSE", {}, ...
    "endErrSSdU", {}, "endErrSSE", {});

for r = 1:size(runDirs, 1)
    runPath = runDirs{r, 1};
    runLabel = runDirs{r, 2};
    files = dir(fullfile(runPath, "out_*.mat"));
    files = sort_struct_by_name(files);
    if isempty(files)
        warning("No files found in %s", runPath);
        continue
    end

    if ~isempty(maxFilesPerRun)
        files = files(1:min(maxFilesPerRun, numel(files)));
    end

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);
        S = load(fpath, "out");
        if ~isfield(S, "out")
            warning("Missing 'out' in %s", fpath);
            continue
        end
        out = S.out;
        if ~isfield(out, "case") || numel(out.case) < 2
            warning("Unexpected case structure in %s", fpath);
            continue
        end

        [ratioSSdU, ratioSSE, tauSSdU, tauSSE, zEval] = compute_ratios(out, cSSdU, cSSE, fullHorizonHours);
        if isempty(ratioSSdU) || isempty(ratioSSE)
            warning("Skipped invalid data in %s", fpath);
            continue
        end

        endRatioSSdU = ratioSSdU(end);
        endRatioSSE = ratioSSE(end);
        endErrSSdU = endRatioSSdU - 1;
        endErrSSE = endRatioSSE - 1;

        records(end+1) = struct( ... %#ok<SAGROW>
            "runLabel", runLabel, ...
            "fileName", files(k).name, ...
            "path", fpath, ...
            "tauSSdU", tauSSdU, ...
            "tauSSE", tauSSE, ...
            "zEval", zEval, ...
            "ratioSSdU", ratioSSdU, ...
            "ratioSSE", ratioSSE, ...
            "endRatioSSdU", endRatioSSdU, ...
            "endRatioSSE", endRatioSSE, ...
            "endErrSSdU", endErrSSdU, ...
            "endErrSSE", endErrSSE);
    end
end

if isempty(records)
    error("No valid run outputs processed.");
end

% ------------------------------------------------------------
% Build matrices for boxplots (one box per out file)
% ------------------------------------------------------------
nRec = numel(records);
dataSSdU = cell(1, nRec);
dataSSE = cell(1, nRec);
zVals = nan(1, nRec);
endErrSSdU = nan(1, nRec);
endErrSSE = nan(1, nRec);

for i = 1:nRec
    dataSSdU{i} = records(i).ratioSSdU(:);
    dataSSE{i} = records(i).ratioSSE(:);
    zVals(i) = records(i).zEval;
    endErrSSdU(i) = records(i).endErrSSdU;
    endErrSSE(i) = records(i).endErrSSE;
end

% ------------------------------------------------------------
% Figures
% ------------------------------------------------------------
fig1 = figure;
subplot(2,1,1);
boxplot(cell2mat(dataSSdU'), group_index(dataSSdU), "symbol", "");
hold on; yline(1, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
ylabel("$\widehat{J}_{TV}(t)/J_{TV,\mathrm{final}}$");
title("\textbf{a}");
grid off; box off
set_font_size(fontSize); format_tick(2, 2);

subplot(2,1,2);
boxplot(cell2mat(dataSSE'), group_index(dataSSE), "symbol", "");
hold on; yline(1, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
ylabel("$\widehat{J}_{track}(t)/J_{track,\mathrm{final}}$");
xlabel("Output file index");
title("\textbf{b}");
grid off; box off
set_font_size(fontSize); format_tick(2, 2);
set_fig_size(1200, 850);

fig2 = figure;
subplot(2,1,1); hold on
for i = 1:nRec
    plot(records(i).tauSSdU, records(i).ratioSSdU, "-", "Color", [plotColors(1,:), 0.25], "LineWidth", 1.0);
end
plot_median_curve(records, "tauSSdU", "ratioSSdU", plotColors(1,:));
yline(1, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("Normalized time $\tau=t/T_f$");
ylabel("$\widehat{J}_{TV}(t)/J_{TV,\mathrm{final}}$");
title("\textbf{a}");
grid off; box off
set_font_size(fontSize); format_tick(2, 2);

subplot(2,1,2); hold on
for i = 1:nRec
    plot(records(i).tauSSE, records(i).ratioSSE, "-", "Color", [plotColors(2,:), 0.25], "LineWidth", 1.0);
end
plot_median_curve(records, "tauSSE", "ratioSSE", plotColors(2,:));
yline(1, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("Normalized time $\tau=t/T_f$");
ylabel("$\widehat{J}_{track}(t)/J_{track,\mathrm{final}}$");
title("\textbf{b}");
grid off; box off
set_font_size(fontSize); format_tick(2, 2);
set_fig_size(1200, 900);

fig3 = figure;
subplot(1,2,1);
scatter(zVals, endErrSSdU, 40, plotColors(1,:), "filled"); hold on
yline(0, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("Evaluation fidelity $z$");
ylabel("End bias $\widehat{J}_{TV}(T_f)/J_{TV,\mathrm{final}} - 1$");
title("\textbf{a}");
grid off; box off
set_font_size(fontSize); format_tick(2, 3);

subplot(1,2,2);
scatter(zVals, endErrSSE, 40, plotColors(2,:), "filled"); hold on
yline(0, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("Evaluation fidelity $z$");
ylabel("End bias $\widehat{J}_{track}(T_f)/J_{track,\mathrm{final}} - 1$");
title("\textbf{b}");
grid off; box off
set_font_size(fontSize); format_tick(2, 3);
set_fig_size(1250, 480);

fig4 = figure;
subplot(1,2,1);
histogram(endErrSSdU, "FaceColor", plotColors(1,:), "FaceAlpha", 0.7, "EdgeColor", "none"); hold on
xline(0, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("End bias in $J_{TV}$ ratio");
ylabel("Count");
title("\textbf{a}");
grid off; box off
set_font_size(fontSize); format_tick(3, 0);

subplot(1,2,2);
histogram(endErrSSE, "FaceColor", plotColors(2,:), "FaceAlpha", 0.7, "EdgeColor", "none"); hold on
xline(0, "--", "Color", plotColors(4,:), "LineWidth", 1.5);
xlabel("End bias in $J_{track}$ ratio");
ylabel("Count");
title("\textbf{b}");
grid off; box off
set_font_size(fontSize); format_tick(3, 0);
set_fig_size(1250, 480);

% ------------------------------------------------------------
% Save outputs
% ------------------------------------------------------------
plotDir = fullfile(projectRoot, "results", "graphical_results");
if ~isfolder(plotDir), mkdir(plotDir); end

print(fig1, fullfile(plotDir, "surrogate_test_boxplots.png"), "-dpng", "-r300");
print(fig1, fullfile(plotDir, "surrogate_test_boxplots.pdf"), "-dpdf", "-bestfit");
print(fig2, fullfile(plotDir, "surrogate_test_ratio_vs_time.png"), "-dpng", "-r300");
print(fig2, fullfile(plotDir, "surrogate_test_ratio_vs_time.pdf"), "-dpdf", "-bestfit");
print(fig3, fullfile(plotDir, "surrogate_test_end_bias_vs_z.png"), "-dpng", "-r300");
print(fig3, fullfile(plotDir, "surrogate_test_end_bias_vs_z.pdf"), "-dpdf", "-bestfit");
print(fig4, fullfile(plotDir, "surrogate_test_end_bias_hist.png"), "-dpng", "-r300");
print(fig4, fullfile(plotDir, "surrogate_test_end_bias_hist.pdf"), "-dpdf", "-bestfit");

% Numeric summary
numDir = fullfile(projectRoot, "results", "numerical results");
if ~isfolder(numDir), mkdir(numDir); end

summaryTxt = fullfile(numDir, "surrogate_test_summary.txt");
fid = fopen(summaryTxt, "w");
if fid == -1
    error("Could not open %s for writing.", summaryTxt);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, "=== Surrogate Extrapolation Validation Summary ===\n");
fprintf(fid, "Processed files: %d\n", nRec);
fprintf(fid, "Full horizon used for fidelity mapping: %.3f h\n\n", fullHorizonHours);

[m1, s1] = mean_std(endErrSSdU);
[m2, s2] = mean_std(endErrSSE);
fprintf(fid, "End-bias stats (ratio-1):\n");
fprintf(fid, "SSdU mean = %.6e, std = %.6e\n", m1, s1);
fprintf(fid, "SSE  mean = %.6e, std = %.6e\n\n", m2, s2);

fprintf(fid, "Per-file end ratios:\n");
fprintf(fid, "run,file,z_eval,end_ratio_SSdU,end_ratio_SSE\n");
for i = 1:nRec
    fprintf(fid, "%s,%s,%.6f,%.8f,%.8f\n", records(i).runLabel, records(i).fileName, ...
        records(i).zEval, records(i).endRatioSSdU, records(i).endRatioSSE);
end

fprintf("Saved surrogate test plots to: %s\n", plotDir);
fprintf("Saved summary to: %s\n", summaryTxt);

% ------------------------------------------------------------
% Helpers
% ------------------------------------------------------------
function [ratioSSdU, ratioSSE, tauSSdU, tauSSE, zEval] = compute_ratios(out, cSSdU, cSSE, fullHorizonHours)
    ratioSSdU = [];
    ratioSSE = [];
    tauSSdU = [];
    tauSSE = [];
    zEval = NaN;

    if ~isfield(out, "case") || numel(out.case) < 2 || ~isfield(out, "T")
        return
    end

    p1_dU = as_col(out.case(1).partial_SSdU);
    p2_dU = as_col(out.case(2).partial_SSdU);
    p1_e = as_col(out.case(1).partial_SSE);
    p2_e = as_col(out.case(2).partial_SSE);
    t = as_col(out.T);

    cumSSdU = cumsum(p1_dU) + cumsum(p2_dU);
    cumSSE = cumsum(p1_e) + cumsum(p2_e);

    fSSdU = t(2:end) / fullHorizonHours;
    fSSE = t / fullHorizonHours;

    tauSSdU = fSSdU;
    tauSSE = fSSE;

    fracSSdU = min(Cheb5(2*fSSdU - 1, cSSdU), 1);
    fracSSE = min(Cheb5(2*fSSE - 1, cSSE), 1);
    fracSSdU = max(fracSSdU, 0.01);
    fracSSE = max(fracSSE, 0.01);

    projSSdU = cumSSdU ./ fracSSdU;
    projSSE = cumSSE ./ fracSSE;

    % Prefer stored final values; fallback to extrapolated final.
    finalSSdU = safe_scalar(out, "SSdU", projSSdU(end));
    finalSSE = safe_scalar(out, "SSE", projSSE(end));

    if ~isfinite(finalSSdU) || finalSSdU <= 0 || ~isfinite(finalSSE) || finalSSE <= 0
        return
    end

    ratioSSdU = projSSdU / finalSSdU;
    ratioSSE = projSSE / finalSSE;

    zEval = out.theta(1);
end

function x = as_col(x)
    x = double(x(:));
end

function y = Cheb5(x, c)
    c = c(:);
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;
    T4 = 8*x.^4 - 8*x.^2 + 1;
    T5 = 16*x.^5 - 20*x.^3 + 5*x;
    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3 + c(5)*T4 + c(6)*T5;
end

function out = group_index(cellData)
    n = numel(cellData);
    out = [];
    for i = 1:n
        out = [out; i*ones(numel(cellData{i}), 1)]; %#ok<AGROW>
    end
end

function plot_median_curve(records, xField, yField, color)
    fGrid = linspace(0, 1, 120)';
    Y = nan(numel(fGrid), numel(records));
    for i = 1:numel(records)
        x = records(i).(xField)(:);
        y = records(i).(yField)(:);
        [xU, idx] = unique(x);
        yU = y(idx);
        if numel(xU) >= 2
            Y(:, i) = interp1(xU, yU, fGrid, "linear", "extrap");
        end
    end
    yMed = median(Y, 2, "omitnan");
    plot(fGrid, yMed, "-", "Color", color, "LineWidth", 2.6);
end

function [ok, cSSdU, cSSE] = try_parse_coeff_txt(pathTxt)
    ok = false;
    cSSdU = [];
    cSSE = [];
    if ~isfile(pathTxt)
        return
    end

    lines = splitlines(string(fileread(pathTxt)));
    idxSSdU = find(contains(lines, "SSdU c0..c5:"), 1, "first");
    idxSSE = find(contains(lines, "SSE c0..c5:"), 1, "first");
    if isempty(idxSSdU) || isempty(idxSSE)
        return
    end

    valsSSdU = nan(6,1);
    valsSSE = nan(6,1);
    for i = 1:6
        valsSSdU(i) = str2double(strtrim(lines(idxSSdU + i)));
        valsSSE(i) = str2double(strtrim(lines(idxSSE + i)));
    end
    if any(~isfinite(valsSSdU)) || any(~isfinite(valsSSE))
        return
    end

    cSSdU = valsSSdU;
    cSSE = valsSSE;
    ok = true;
end

function filesOut = sort_struct_by_name(filesIn)
    [~, idx] = sort({filesIn.name});
    filesOut = filesIn(idx);
end

function v = safe_scalar(S, fieldName, fallback)
    if isfield(S, fieldName)
        v = double(S.(fieldName));
        if numel(v) ~= 1 || ~isfinite(v)
            v = fallback;
        end
    else
        v = fallback;
    end
end

function [m, s] = mean_std(x)
    x = x(isfinite(x));
    if isempty(x)
        m = NaN; s = NaN;
    else
        m = mean(x);
        s = std(x);
    end
end
