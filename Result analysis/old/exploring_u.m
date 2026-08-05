%% Explore reduced-variation input profiles for the selected controller
% For node points spaced every ceil(0.1/dt) samples, enforce
%   sum(U_star(1:k,:)) = sum(U(1:k,:))
% at each node k and minimize the total variation of U_star.
%
% Here the solve is done explicitly, block by block. For each interval
% between consecutive nodes and for each input channel, solve the LP
%   min sum_i |u_i - u_{i-1}|
% subject to
%   sum_i u_i = sum_i u_i(original block),
% with the first increment measured against the previous optimized sample.

clear; close all; clc

projectRoot = fileparts(fileparts(mfilename("fullpath")));
resultsRoot = fullfile(projectRoot, "results");

cfg = struct();
cfg.schedule_root = fullfile(resultsRoot, "setpoint_schedule_xsp_7_13_16");
cfg.scenario = "same_noise"; % "same_noise" or "no_noise"
cfg.selected_controller_id = "ts_20260211_122653_modified";
cfg.node_spacing_h = 0.1;
cfg.graphics_dir = fullfile(resultsRoot, "graphics");
cfg.save_outputs = true;
cfg.input_lb = 0;
cfg.input_ub = 0.4;

if ~isfolder(cfg.graphics_dir)
    mkdir(cfg.graphics_dir);
end

runDir = fullfile(cfg.schedule_root, cfg.scenario);
matPath = fullfile(runDir, "out_schedule_" + cfg.selected_controller_id + ".mat");
if ~isfile(matPath)
    error("Selected-controller schedule file not found: %s", matPath);
end

S = load(matPath, "out", "ctrl");
if ~isfield(S, "out") || ~isfield(S.out, "case") || isempty(S.out.case)
    error("Missing out.case in %s", matPath);
end

out = S.out;
nCase = min(numel(out.case), 2);
plotColors = nature_methods_colors(3);
inputLabels = ["$F_{in}$ [L/h]", "$F_m$ [L/h]", "$F_{out}$ [L/h]"];
casePanelLabels = ["a", "b"];
fontSize = 14;

summaryRows = table();
solutions = struct("t", {}, "U", {}, "Ustar", {}, "node_idx", {}, "node_t", {}, "dt", {});

for c = 1:nCase
    caseData = out.case(c);
    U = get_case_inputs(caseData);
    dt = get_case_dt(caseData, size(U, 1));
    t = (0:size(U, 1)-1).' * dt;

    nodeStride = max(1, ceil(cfg.node_spacing_h / dt));
    [Ustar, nodeIdx] = solve_blockwise_tv_linprog(U, nodeStride, cfg.input_lb, cfg.input_ub);
    nodeTime = t(nodeIdx);

    cumU = cumsum(U, 1);
    cumUstar = cumsum(Ustar, 1);
    nodeResidual = cumUstar(nodeIdx, :) - cumU(nodeIdx, :);
    tvU = compute_total_variation(U);
    tvUstar = compute_total_variation(Ustar);
    maxDiffU = compute_max_abs_diff(U);
    maxDiffUstar = compute_max_abs_diff(Ustar);

    summaryRows = [summaryRows; table( ...
        c, dt, nodeStride, numel(nodeIdx), max(abs(nodeResidual), [], "all"), ...
        tvU(1), tvU(2), tvU(3), tvUstar(1), tvUstar(2), tvUstar(3), ...
        maxDiffU(1), maxDiffU(2), maxDiffU(3), ...
        maxDiffUstar(1), maxDiffUstar(2), maxDiffUstar(3), ...
        'VariableNames', { ...
        'case_id', 'dt_h', 'node_stride', 'n_nodes', 'max_node_sum_error', ...
        'TV_Fin', 'TV_Fm', 'TV_Fout', 'TVstar_Fin', 'TVstar_Fm', 'TVstar_Fout', ...
        'maxDiff_Fin', 'maxDiff_Fm', 'maxDiff_Fout', ...
        'maxDiffStar_Fin', 'maxDiffStar_Fm', 'maxDiffStar_Fout'})]; %#ok<AGROW>

    solutions(c).t = t;
    solutions(c).U = U;
    solutions(c).Ustar = Ustar;
    solutions(c).node_idx = nodeIdx;
    solutions(c).node_t = nodeTime;
    solutions(c).dt = dt;
end

disp("Reduced-TV input reconstruction summary:")
disp(summaryRows)

fig = figure("Color", "w", "Name", "Selected controller: U vs U*");
tiledlayout(fig, nCase, 3, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [120 120 1400 760]);

for c = 1:nCase
    t = solutions(c).t;
    U = solutions(c).U;
    Ustar = solutions(c).Ustar;

    for u = 1:min(3, size(U, 2))
        ax = nexttile((c - 1) * 3 + u); hold(ax, "on");
        stairs(ax, t, U(:, u), "-", "Color", [0 0 0], "LineWidth", 1.0, "DisplayName", "$U$");
        stairs(ax, t, Ustar(:, u), "--", "Color", plotColors(1, :), "LineWidth", 2.0, "DisplayName", "$U^\star$");

        if c < nCase
            xlabelText = "";
        else
            xlabelText = "Time [h]";
        end

        xlabel(ax, xlabelText, "Interpreter", "latex");
        ylabel(ax, inputLabels(u), "Interpreter", "latex");
        if u == 1
            title(ax, sprintf('$\\mathbf{%s}$', char(casePanelLabels(c))), "Interpreter", "latex");
            ax.TitleHorizontalAlignment = "left";
        else
            title(ax, "");
        end
        if c == 1 && u == 1
            legend(ax, "Interpreter", "latex", "Location", "best", "Box", "off");
        end
        ax.TickLabelInterpreter = "latex";
        set(ax, "FontSize", fontSize);
        grid(ax, "off");
        box(ax, "off");
    end
end

if cfg.save_outputs
    outStem = fullfile(cfg.graphics_dir, "exploring_u_selected_controller_" + cfg.scenario);
    exportgraphics(fig, outStem + ".png", "Resolution", 300);
    exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
    fprintf("Saved: %s\n", outStem + ".png");
    fprintf("Saved: %s\n", outStem + ".pdf");
end


function U = get_case_inputs(caseData)
if ~isfield(caseData, "U") || isempty(caseData.U)
    error("Missing caseData.U in selected-controller schedule results.");
end
U = double(caseData.U);
if size(U, 2) < 3
    error("Expected three input channels, found %d.", size(U, 2));
end
U = U(:, 1:3);
end

function dt = get_case_dt(caseData, nSample)
if isfield(caseData, "dt")
    dt = double(caseData.dt);
elseif isfield(caseData, "tf")
    dt = double(caseData.tf) / max(nSample - 1, 1);
else
    dt = 1/60;
end
end

function [Ustar, nodeIdx] = solve_blockwise_tv_linprog(U, nodeStride, inputLb, inputUb)
[nSample, nInput] = size(U);
nodeIdx = unique([nodeStride:nodeStride:nSample, nSample]);
blockStart = [1, nodeIdx(1:end-1) + 1];
blockEnd = nodeIdx;

Ustar = zeros(nSample, nInput);
linprogOptions = optimoptions("linprog", "Display", "none", "Algorithm", "dual-simplex");

for b = 1:numel(nodeIdx)
    idx = blockStart(b):blockEnd(b);
    for j = 1:nInput
        if b == 1
            uPrev = U(1, j);
        else
            uPrev = Ustar(blockStart(b) - 1, j);
        end
        Ustar(idx, j) = solve_single_block_lp(U(idx, j), uPrev, inputLb, inputUb, linprogOptions);
    end
end
end

function uBlockStar = solve_single_block_lp(uBlock, uPrev, inputLb, inputUb, linprogOptions)
m = numel(uBlock);
sumTarget = sum(uBlock);

if sumTarget < m * inputLb - 1e-10 || sumTarget > m * inputUb + 1e-10
    error("Block LP infeasible under bounds: length=%d, target sum=%.6g, feasible range=[%.6g, %.6g].", ...
        m, sumTarget, m * inputLb, m * inputUb);
end

if m == 1
    uBlockStar = sumTarget;
    return
end

nVar = 2 * m;
f = [zeros(m, 1); ones(m, 1)];

Aeq = zeros(1, nVar);
Aeq(1, 1:m) = 1;
beq = sumTarget;

A = zeros(2 * m, nVar);
b = zeros(2 * m, 1);

% First increment is measured relative to the previous optimized sample.
A(1, 1) = 1;
A(1, m + 1) = -1;
b(1) = uPrev;

A(2, 1) = -1;
A(2, m + 1) = -1;
b(2) = -uPrev;

for i = 2:m
    rowPos = 2 * i - 1;
    rowNeg = 2 * i;

    A(rowPos, i) = 1;
    A(rowPos, i - 1) = -1;
    A(rowPos, m + i) = -1;

    A(rowNeg, i) = -1;
    A(rowNeg, i - 1) = 1;
    A(rowNeg, m + i) = -1;
end

lb = [inputLb * ones(m, 1); zeros(m, 1)];
ub = [inputUb * ones(m, 1); inf(m, 1)];

[xOpt, ~, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub, linprogOptions);
if exitflag <= 0 || isempty(xOpt)
    error("Block LP failed: exitflag=%d, block length=%d, target sum=%.6g.", exitflag, m, sumTarget);
end

uBlockStar = xOpt(1:m);
end

function tvVal = compute_total_variation(U)
if size(U, 1) <= 1
    tvVal = zeros(1, size(U, 2));
    return
end
tvVal = sum(abs(diff(U, 1, 1)), 1);
end

function maxDiffVal = compute_max_abs_diff(U)
if size(U, 1) <= 1
    maxDiffVal = zeros(1, size(U, 2));
    return
end
maxDiffVal = max(abs(diff(U, 1, 1)), [], 1);
end
