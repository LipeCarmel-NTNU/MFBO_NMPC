function plot_refined_frontier_change()
%PLOT_REFINED_FRONTIER_CHANGE Visualize how the Pareto frontier changed after refinement.
%
% Inputs (from disk):
%   - Original runs:
%       results/run1/results.csv
%       results/run2/results.csv
%   - Refined full-fidelity evaluations:
%       results/final_fidelity_same_noise/run1_full_f1_same_noise/results_full.csv
%       results/final_fidelity_same_noise/run2_full_f1_same_noise/results_full.csv
%
% Outputs:
%   - results/graphical_results/refined_frontier_change.png
%   - results/graphical_results/refined_frontier_change.pdf
%
% Figure content:
%   - Original non-Pareto samples only (gray), excluding DOE points
%     (iterations 1..20 in run1/run2)
%   - Refined points (blue)
%   - Refined Pareto frontier (green)
%   - Promoted points with original z < 1 (red stars)
%   - Thin connector lines for promoted points (original -> refined)

    close all; clc
    scriptDir = fileparts(mfilename("fullpath"));
    projectRoot = fileparts(scriptDir);
    addpath(genpath(fullfile(projectRoot, "dependencies")));

    originalFiles = [
        struct("runLabel","Case 1", "path", fullfile(projectRoot, "results", "run1", "results.csv"));
        struct("runLabel","Case 2", "path", fullfile(projectRoot, "results", "run2", "results.csv"))
    ];

    refinedFiles = [
        struct("runLabel","Case 1 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run1_full_f1_same_noise", "results_full.csv"));
        struct("runLabel","Case 2 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run2_full_f1_same_noise", "results_full.csv"))
    ];

    doeIterationsPerRun = 20;
    print_runtime_cfg(originalFiles, refinedFiles, doeIterationsPerRun);

    T_orig = load_and_stack_csv(originalFiles, doeIterationsPerRun);
    T_ref = load_and_stack_csv(refinedFiles, 0);

    if isempty(T_orig)
        error("No rows loaded from original result files.");
    end
    if isempty(T_ref)
        error("No rows loaded from refined result files.");
    end

    T_jtv = compute_jtv_change_table(T_orig, T_ref);
    fprintf("\n=== J_TV change (original -> refined), sorted by |delta| descending ===\n");
    disp(T_jtv);
    print_top1_cfg(projectRoot, T_jtv);

    % Combined Pareto masks.
    isParetoOrig = compute_pareto_mask(double(T_orig.SSE), double(T_orig.SSdU));
    isParetoRef = compute_pareto_mask(double(T_ref.SSE), double(T_ref.SSdU));

    % Timestamp-level matching between original and refined points.
    [commonTs, idxOrig, idxRef] = intersect(string(T_orig.timestamp), string(T_ref.timestamp), "stable"); %#ok<ASGLU>

    C = nature_methods_colors();
    colOrigAll = 0.80 * [1 1 1];
    colRefAll = C.SkyBlue;
    colRefPareto = C.BluishGreen;
    colShift = 0.35 * [1 1 1];
    colPromoted = C.Vermillion;

    % Promotion definition (discard original Pareto points):
    % promoted = not Pareto in original AND Pareto in refined.
    origParetoMatched = isParetoOrig(idxOrig);
    refParetoMatched = isParetoRef(idxRef);
    zOrigMatched = double(T_orig.z_eval(idxOrig));
    promotedMask = ~origParetoMatched & refParetoMatched;
    promotedZlt1Mask = promotedMask & isfinite(zOrigMatched) & (zOrigMatched < 1 - 1e-12);
    promotedIdxOrig = idxOrig(promotedZlt1Mask);
    promotedIdxRef = idxRef(promotedZlt1Mask);
    promotedTs = string(T_orig.timestamp(promotedIdxOrig));

    fig = figure("Color", "w", "Toolbar", "none");
    ax = axes(fig); hold(ax, "on");

    % Original non-Pareto samples only (original Pareto discarded).
    scatter(ax, double(T_orig.SSdU(~isParetoOrig)), double(T_orig.SSE(~isParetoOrig)), 24, ...
        "filled", "MarkerFaceColor", colOrigAll, "MarkerEdgeColor", "none", ...
        "DisplayName", "Original non-Pareto");

    % Refined samples and refined Pareto.
    scatter(ax, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
        "filled", "MarkerFaceColor", colRefAll, "MarkerEdgeColor", "none", ...
        "DisplayName", "Refined samples");

    draw_frontier(ax, T_ref(isParetoRef, :), colRefPareto, "Refined Pareto");

    % Connect only promoted points (original non-Pareto -> refined Pareto with z<1).
    for i = 1:numel(promotedIdxOrig)
        x = [double(T_orig.SSdU(promotedIdxOrig(i))), double(T_ref.SSdU(promotedIdxRef(i)))];
        y = [double(T_orig.SSE(promotedIdxOrig(i))), double(T_ref.SSE(promotedIdxRef(i)))];
        plot(ax, x, y, "-", "Color", colShift, "LineWidth", 1.0, "HandleVisibility", "off");
    end

    if ~isempty(promotedIdxRef)
        scatter(ax, double(T_ref.SSdU(promotedIdxRef)), double(T_ref.SSE(promotedIdxRef)), 140, ...
            "p", "MarkerFaceColor", colPromoted, "MarkerEdgeColor", "k", "LineWidth", 0.7, ...
            "DisplayName", "Promoted (z < 1)");
    end

    set(ax, "XScale", "log", "YScale", "log");
    xlabel(ax, "$J_{\mathrm{TV}}$", "Interpreter", "latex");
    ylabel(ax, "$J_{\mathrm{track}}$", "Interpreter", "latex");
    title(ax, "$\mathbf{Pareto~Frontier~Change:~Original~vs~Refined}$", "Interpreter", "latex");
    grid(ax, "off");
    box(ax, "off");
    set(ax, "FontSize", 16);
    format_tick(1, 1);

    lg = legend(ax, "Location", "southwest");
    lg.Interpreter = "latex";
    lg.Box = "off";

    txt = sprintf(['matched timestamps: %d\n' ...
        'orig Pareto points: %d\n' ...
        'refined Pareto points: %d\n' ...
        'promoted to Pareto (z < 1): %d'], ...
        numel(idxOrig), nnz(isParetoOrig), nnz(isParetoRef), numel(promotedIdxRef));
    text(ax, 0.98, 0.04, txt, "Units", "normalized", ...
        "HorizontalAlignment", "right", "VerticalAlignment", "bottom", ...
        "Interpreter", "none", "FontSize", 12);

    outDir = fullfile(projectRoot, "results", "graphical_results");
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    outStem = fullfile(outDir, "refined_frontier_change");
    exportgraphics(fig, outStem + ".png", "Resolution", 300);
    exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");

    outTxtDir = fullfile(projectRoot, "results", "txt results");
    if ~isfolder(outTxtDir)
        mkdir(outTxtDir);
    end
    reportPath = fullfile(outTxtDir, "refined_promoted_frontier_z_lt_1.txt");
    write_promoted_report(reportPath, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef);

    fprintf("Saved: %s\n", outStem + ".png");
    fprintf("Saved: %s\n", outStem + ".pdf");
    fprintf("Saved: %s\n", reportPath);
    fprintf("J_TV rows compared: %d\n", height(T_jtv));
    fprintf("Matched timestamps: %d\n", numel(idxOrig));
    fprintf("Original Pareto points: %d\n", nnz(isParetoOrig));
    fprintf("Refined Pareto points: %d\n", nnz(isParetoRef));
    fprintf("Promoted to Pareto with z < 1: %d\n", numel(promotedIdxRef));
end


function print_runtime_cfg(originalFiles, refinedFiles, doeIterationsPerRun)
%PRINT_RUNTIME_CFG Print analysis/runtime settings to command window.
    fprintf("=== plot_refined_frontier_change settings ===\n");
    fprintf("MATLAB version: %s\n", version);
    fprintf("Current folder: %s\n", pwd);
    fprintf("DOE exclusion in original runs: drop first %d iterations per run\n", doeIterationsPerRun);

    p = gcp("nocreate");
    if isempty(p)
        fprintf("Parallel pool: not running\n");
        workers = 0;
    else
        workers = p.NumWorkers;
        fprintf("Parallel pool: running (%d workers)\n", workers);
    end

    fprintf("Configured worker count (observed): %d\n", workers);
    fprintf("Original CSV sources:\n");
    for k = 1:numel(originalFiles)
        fprintf("  - %s\n", originalFiles(k).path);
    end
    fprintf("Refined CSV sources:\n");
    for k = 1:numel(refinedFiles)
        fprintf("  - %s\n", refinedFiles(k).path);
    end
end


function T = compute_jtv_change_table(T_orig, T_ref)
%COMPUTE_JTV_CHANGE_TABLE Build matched timestamp table and sort by |delta J_TV|.
    [commonTs, idxOrig, idxRef] = intersect(string(T_orig.timestamp), string(T_ref.timestamp), "stable");
    if isempty(commonTs)
        error("No matching timestamps between original and refined tables.");
    end
    T = table(commonTs, ...
        double(T_orig.SSdU(idxOrig)), ...
        double(T_ref.SSdU(idxRef)), ...
        'VariableNames', {'timestamp', 'JTV_original', 'JTV_refined'});
    T.delta_JTV = T.JTV_refined - T.JTV_original;
    T.delta_pct = 100 * T.delta_JTV ./ T.JTV_original;
    T.abs_delta_JTV = abs(T.delta_JTV);
    T = sortrows(T, "abs_delta_JTV", "descend");
end


function print_top1_cfg(projectRoot, T_jtv)
%PRINT_TOP1_CFG Print original/refined cfg for top-1 |delta J_TV| timestamp.
    if isempty(T_jtv)
        fprintf("\nNo J_TV rows available for cfg comparison.\n");
        return
    end
    ts = string(T_jtv.timestamp(1));

    origMat = find_mat_for_timestamp(projectRoot, ts, true);
    refMat = find_mat_for_timestamp(projectRoot, ts, false);
    if strlength(origMat) == 0 || strlength(refMat) == 0
        fprintf("\nTop-1 cfg print skipped (missing MAT files) for timestamp %s.\n", ts);
        return
    end

    Sorig = load(origMat);
    Sref = load(refMat);

    fprintf("\n=== Top-1 timestamp by |delta J_TV|: %s ===\n", ts);
    fprintf("\n--- Original cfg (%s) ---\n", origMat);
    print_cfg_struct(Sorig);
    disp(Sorig.cfg_run)
    disp(Sorig.base)

    fprintf("\n--- Refined cfg (%s) ---\n", refMat);
    print_cfg_struct(Sref);
    disp(Sref.cfg_run)
    disp(Sref.base)
end


function matPath = find_mat_for_timestamp(projectRoot, ts, isOriginal)
%FIND_MAT_FOR_TIMESTAMP Locate original/refined MAT file for a timestamp.
    if isOriginal
        cands = [
            fullfile(projectRoot, "results", "run1", "out_" + ts + ".mat");
            fullfile(projectRoot, "results", "run2", "out_" + ts + ".mat")
        ];
    else
        cands = [
            fullfile(projectRoot, "results", "final_fidelity_same_noise", "run1_full_f1_same_noise", "out_full_" + ts + ".mat");
            fullfile(projectRoot, "results", "final_fidelity_same_noise", "run2_full_f1_same_noise", "out_full_" + ts + ".mat")
        ];
    end
    matPath = "";
    for i = 1:numel(cands)
        if isfile(cands(i))
            matPath = string(cands(i));
            return
        end
    end
end


function print_cfg_struct(S)
%PRINT_CFG_STRUCT Display cfg/theta/sigma_y values if present.
    if isfield(S, "out") && isfield(S.out, "cfg")
        cfg = S.out.cfg;
        if isfield(cfg, "f"), fprintf("f = %.12g\n", cfg.f); end
        if isfield(cfg, "m"), fprintf("m = %d\n", cfg.m); end
        if isfield(cfg, "p"), fprintf("p = %d\n", cfg.p); end
        if isfield(cfg, "Q"), fprintf("Q =\n"); disp(cfg.Q); end
        if isfield(cfg, "Ru"), fprintf("Ru =\n"); disp(cfg.Ru); end
        if isfield(cfg, "Rdu"), fprintf("Rdu =\n"); disp(cfg.Rdu); end
    else
        fprintf("No out.cfg found.\n");
    end
    if isfield(S, "out") && isfield(S.out, "theta")
        fprintf("theta =\n");
        disp(S.out.theta);
    end
    if isfield(S, "cfg_run") && isfield(S.cfg_run, "sigma_y")
        fprintf("sigma_y = ");
        disp(S.cfg_run.sigma_y);
    end
end


function T = load_and_stack_csv(fileDefs, dropFirstNRows)
%LOAD_AND_STACK_CSV Load multiple CSV files and stack relevant columns.
    if nargin < 2
        dropFirstNRows = 0;
    end
    T = table();
    for k = 1:numel(fileDefs)
        p = fileDefs(k).path;
        if ~isfile(p)
            warning("Missing file: %s", p);
            continue
        end
        Tk = readtable(p, "TextType", "string");
        required = ["timestamp", "SSE", "SSdU"];
        for c = required
            if ~ismember(c, string(Tk.Properties.VariableNames))
                error("Missing column '%s' in %s", c, p);
            end
        end
        Tk.run_label = repmat(string(fileDefs(k).runLabel), height(Tk), 1);
        Tk.iteration = (1:height(Tk)).';
        if dropFirstNRows > 0
            Tk = Tk(Tk.iteration > dropFirstNRows, :);
        end
        if ismember("theta_1", string(Tk.Properties.VariableNames))
            Tk.z_eval = double(Tk.theta_1);
        else
            Tk.z_eval = nan(height(Tk), 1);
        end
        Tk = Tk(:, ["run_label", "iteration", "timestamp", "SSE", "SSdU", "z_eval"]);
        T = [T; Tk]; %#ok<AGROW>
    end
end


function draw_frontier(ax, Tpareto, colorVal, legendLabel)
%DRAW_FRONTIER Draw sorted Pareto points and connecting line.
    if isempty(Tpareto)
        return
    end
    x = double(Tpareto.SSdU);
    y = double(Tpareto.SSE);
    [xSort, idx] = sort(x, "ascend");
    ySort = y(idx);

    plot(ax, xSort, ySort, "-", "LineWidth", 2.4, "Color", colorVal, "DisplayName", legendLabel);
    scatter(ax, xSort, ySort, 82, ...
        "o", "MarkerFaceColor", "none", "MarkerEdgeColor", colorVal, ...
        "LineWidth", 1.8, "HandleVisibility", "off");
end


function isPareto = compute_pareto_mask(J_track, J_TV)
%COMPUTE_PARETO_MASK Return non-dominated mask for minimization objectives.
    n = numel(J_track);
    isPareto = true(n, 1);
    for i = 1:n
        dominated = (J_track <= J_track(i)) & (J_TV <= J_TV(i)) & ...
                    ((J_track < J_track(i)) | (J_TV < J_TV(i)));
        dominated(i) = false;
        if any(dominated)
            isPareto(i) = false;
        end
    end
end


function write_promoted_report(reportPath, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef)
%WRITE_PROMOTED_REPORT Save promoted refined-Pareto points with original z<1.
    fid = fopen(reportPath, "w");
    if fid < 0
        warning("Could not write promoted-point report: %s", reportPath);
        return
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "Promoted to refined Pareto frontier with original z < 1\n");
    fprintf(fid, "Original dataset excludes DOE rows (iterations 1..20 in run1/run2).\n");
    fprintf(fid, "Definition: original non-Pareto -> refined Pareto AND original theta_1 < 1\n");
    fprintf(fid, "Count: %d\n\n", numel(promotedTs));
    fprintf(fid, "timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n");

    for i = 1:numel(promotedTs)
        io = promotedIdxOrig(i);
        ir = promotedIdxRef(i);
        fprintf(fid, "%s,%.17g,%.17g,%.17g,%.17g,%.17g\n", ...
            promotedTs(i), ...
            double(T_orig.SSE(io)), double(T_orig.SSdU(io)), double(T_orig.z_eval(io)), ...
            double(T_ref.SSE(ir)), double(T_ref.SSdU(ir)));
    end
end
