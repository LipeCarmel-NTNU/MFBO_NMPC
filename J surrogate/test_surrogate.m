% test_surrogate.m
% Validate Chebyshev-based multifidelity extrapolation against stored run outputs.
% Kept output: normalized prediction trajectories only.

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
fullHorizonHours = 10.0;  % T in tau = t / (z*T)
maxFilesPerRun = 10;      % set [] to include all files
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
    "endRatioSSdU", {}, "endRatioSSE", {});

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

        records(end+1) = struct( ... %#ok<SAGROW>
            "runLabel", runLabel, ...
            "fileName", files(k).name, ...
            "path", fpath, ...
            "tauSSdU", tauSSdU, ...
            "tauSSE", tauSSE, ...
            "zEval", zEval, ...
            "ratioSSdU", ratioSSdU, ...
            "ratioSSE", ratioSSE, ...
            "endRatioSSdU", ratioSSdU(end), ...
            "endRatioSSE", ratioSSE(end));
    end
end

if isempty(records)
    error("No valid run outputs processed.");
end

nRec = numel(records);

% ------------------------------------------------------------
% Figure: normalized predicted cost trajectories
% ------------------------------------------------------------
fig = figure;
subplot(2,1,1); hold on
for i = 1:nRec
    plot(records(i).tauSSdU, records(i).ratioSSdU, "-", "Color", [plotColors(1,:), 0.25], "LineWidth", 1.0);
end
plot_median_curve(records, "tauSSdU", "ratioSSdU", plotColors(1,:));
yline(1, "--", "Color", [0.2 0.2 0.2], "LineWidth", 1.8);
xlabel("Normalized time $\tau=t/(z\cdot T)$");
ylabel("$\widehat{J}_{TV}(t)/J_{TV,\mathrm{final}}$");
title("\textbf{a}");
xlim([0, 1]); ylim([0, 1.5]);
grid off; box off
set_font_size(fontSize); format_tick(2, 2);

subplot(2,1,2); hold on
for i = 1:nRec
    plot(records(i).tauSSE, records(i).ratioSSE, "-", "Color", [plotColors(2,:), 0.25], "LineWidth", 1.0);
end
plot_median_curve(records, "tauSSE", "ratioSSE", plotColors(2,:));
yline(1, "--", "Color", [0.2 0.2 0.2], "LineWidth", 1.8);
xlabel("Normalized time $\tau=t/(z\cdot T)$");
ylabel("$\widehat{J}_{track}(t)/J_{track,\mathrm{final}}$");
title("\textbf{b}");
xlim([0, 1]); ylim([0, 1.5]);
grid off; box off
set_font_size(fontSize); format_tick(2, 2);
set_fig_size(1200, 900);

% ------------------------------------------------------------
% Save outputs
% ------------------------------------------------------------
plotDir = fullfile(projectRoot, "results", "graphical_results");
if ~isfolder(plotDir), mkdir(plotDir); end

print(fig, fullfile(plotDir, "surrogate_test_ratio_vs_time.png"), "-dpng", "-r300");
print(fig, fullfile(plotDir, "surrogate_test_ratio_vs_time.pdf"), "-dpdf", "-bestfit");

% Numeric summary
numDir = fullfile(projectRoot, "results", "numerical results");
if ~isfolder(numDir), mkdir(numDir); end

summaryTxt = fullfile(numDir, "surrogate_test_summary.txt");
fid = fopen(summaryTxt, "w");
if fid == -1
    error("Could not open %s for writing.", summaryTxt);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, "=== Surrogate Extrapolation Validation Summary ===\n");
fprintf(fid, "Processed files: %d\n", nRec);
fprintf(fid, "Full horizon used for fidelity mapping: %.3f h\n\n", fullHorizonHours);
fprintf(fid, "Per-file end ratios:\n");
fprintf(fid, "run,file,z_eval,end_ratio_SSdU,end_ratio_SSE\n");
for i = 1:nRec
    fprintf(fid, "%s,%s,%.6f,%.8f,%.8f\n", records(i).runLabel, records(i).fileName, ...
        records(i).zEval, records(i).endRatioSSdU, records(i).endRatioSSE);
end

fprintf("Saved surrogate test plot to: %s\n", plotDir);
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

    % partial_SSdU starts at iteration 2; partial_SSE starts at iteration 1.
    NdU = min([numel(p1_dU), numel(p2_dU), max(numel(t) - 1, 0)]);
    NSE = min([numel(p1_e), numel(p2_e), numel(t)]);
    if NdU < 3 || NSE < 3
        return
    end

    p1_dU = p1_dU(1:NdU); p2_dU = p2_dU(1:NdU);
    p1_e = p1_e(1:NSE); p2_e = p2_e(1:NSE);
    tSSdU = t(2:(NdU + 1));
    tSSE = t(1:NSE);

    cumSSdU = cumsum(p1_dU) + cumsum(p2_dU);
    cumSSE = cumsum(p1_e) + cumsum(p2_e);

    fSSdU = tSSdU / fullHorizonHours;
    fSSE = tSSE / fullHorizonHours;
    if any(fSSdU < 0 | fSSdU > 1) || any(fSSE < 0 | fSSE > 1)
        error("Encountered fidelity values outside [0,1] while building surrogate test ratios.");
    end

    fracSSdU = min(Cheb5(2*fSSdU - 1, cSSdU), 1);
    fracSSE = min(Cheb5(2*fSSE - 1, cSSE), 1);
    fracSSdU = max(fracSSdU, 0.01);
    fracSSE = max(fracSSE, 0.01);

    projSSdU = cumSSdU ./ fracSSdU;
    projSSE = cumSSE ./ fracSSE;

    finalSSdU = safe_scalar(out, "SSdU", projSSdU(end));
    finalSSE = safe_scalar(out, "SSE", projSSE(end));
    if ~isfinite(finalSSdU) || finalSSdU <= 0 || ~isfinite(finalSSE) || finalSSE <= 0
        return
    end

    ratioSSdU = projSSdU / finalSSdU;
    ratioSSE = projSSE / finalSSE;

    tfHours = safe_scalar(out, "tf", tSSE(end));
    zEval = min(max(tfHours / fullHorizonHours, 0), 1);
    tauSSdU = tSSdU / max(tfHours, eps);
    tauSSE = tSSE / max(tfHours, eps);
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
