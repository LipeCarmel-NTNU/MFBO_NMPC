function plot_refined_frontier_change()
%PLOT_REFINED_FRONTIER_CHANGE Compare original vs refined frontier decisions.
%
% Intent
% - Answer one question: "How did the candidate controllers chosen by the
%   optimization frontier change after full-fidelity reevaluation?"
% - Keep the comparison fair by excluding DOE points from the original
%   optimization history.
% - Compare only like-for-like controller identities (run + timestamp).
%
% Data contract
% - Original optimization histories:
%   results/run1/results.csv, results/run2/results.csv
% - Refined reevaluations:
%   results/final_fidelity_same_noise/run1_full_f1_same_noise/results_full.csv
%   results/final_fidelity_same_noise/run2_full_f1_same_noise/results_full.csv
%
% Figure intent
% - Panel a: reference picture of the combined optimization frontier from
%   the matched non-DOE points (context view).
% - Panel b: direct before/after comparison where:
%   - original points remain visible by run (marker + color identity),
%   - refined evaluations are overlaid,
%   - refined Pareto points are highlighted (purple diamonds).
%
% Filtering intent (fairness guards)
% - DOE exclusion: original iterations 1..20 are never eligible.
% - Eligibility: refined points are only kept if they map to an original
%   non-DOE Pareto controller identity.
% - Identity matching: (run_key, timestamp), not timestamp alone.
%
% Outputs
% - results/graphical_results/refined_frontier_change.png/.pdf
% - results/txt results/refined_promoted_frontier_z_lt_1.txt

    close all; clc
    scriptDir = fileparts(mfilename("fullpath"));
    projectRoot = fileparts(scriptDir);
    addpath(genpath(fullfile(projectRoot, "dependencies")));

    originalFiles = [
        struct("runKey","run1", "runLabel","Case 1", "path", fullfile(projectRoot, "results", "run1", "results.csv"));
        struct("runKey","run2", "runLabel","Case 2", "path", fullfile(projectRoot, "results", "run2", "results.csv"))
    ];

    refinedFiles = [
        struct("runKey","run1", "runLabel","Case 1 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run1_full_f1_same_noise", "results_full.csv"));
        struct("runKey","run2", "runLabel","Case 2 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run2_full_f1_same_noise", "results_full.csv"))
    ];

    doeIterationsPerRun = 20;
    print_runtime_cfg(originalFiles, refinedFiles, doeIterationsPerRun);

    T_orig_full = load_and_stack_csv(originalFiles, 0);
    T_orig = T_orig_full(T_orig_full.iteration > doeIterationsPerRun, :);
    T_ref = load_and_stack_csv(refinedFiles, 0);

    if isempty(T_orig_full)
        error("No rows loaded from original result files.");
    end
    if isempty(T_orig)
        error("No original rows remain after DOE filtering.");
    end
    if isempty(T_ref)
        error("No rows loaded from refined result files.");
    end

    % Original Pareto mask on optimization-only set (DOE already removed).
    isParetoOrig = compute_pareto_mask(double(T_orig.SSE), double(T_orig.SSdU));

    % Keep only refined points that map to original non-DOE Pareto points.
    [T_ref, refFilterInfo] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig, T_orig_full);
    print_refined_filter_summary(refFilterInfo);
    if isempty(T_ref)
        error("No refined rows remain after filtering to original non-DOE Pareto points.");
    end

    T_jtv = compute_jtv_change_table(T_orig, T_ref);
    fprintf("\n=== J_TV change (original -> refined), sorted by |delta| descending ===\n");
    disp(T_jtv);
    print_top1_cfg(projectRoot, T_jtv);

    % Combined Pareto masks.
    isParetoRef = compute_pareto_mask(double(T_ref.SSE), double(T_ref.SSdU));

    % Timestamp-level matching between original and refined points.
    [commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), "stable"); %#ok<ASGLU>
    [matchedRunKey, matchedTs] = split_match_key(commonKeys);

    NATURE_COLOR = nature_methods_colors();
    plotColors = nature_methods_colors(3); % Blue, BluishGreen, ReddishPurple
    colRefAll = [0.60 0.82 0.98];
    colPromoted = plotColors(2, :);

    % Left-panel data: all combined optimization samples + matched subset frontier.
    T1 = T_orig(T_orig.run_key == "run1", :);
    T2 = T_orig(T_orig.run_key == "run2", :);
    Tall = [T1; T2];
    matchedKeySet = string(commonKeys);
    isLeft1 = ismember(string(T1.match_key), matchedKeySet);
    isLeft2 = ismember(string(T2.match_key), matchedKeySet);
    TleftFrontPool = [T1(isLeft1, :); T2(isLeft2, :)];
    leftMask = compute_pareto_mask(double(TleftFrontPool.SSE), double(TleftFrontPool.SSdU));
    Tf = TleftFrontPool(leftMask, :);

    % Safety check: left frontier must not include DOE-derived points.
    doeKeys = string(T_orig_full.match_key(T_orig_full.iteration <= doeIterationsPerRun));
    if any(ismember(string(Tf.match_key), doeKeys))
        error("Left-panel frontier unexpectedly contains DOE-derived points.");
    end

    % Promotion definition (discard original Pareto points):
    % promoted = not Pareto in original AND Pareto in refined.
    origParetoMatched = isParetoOrig(idxOrig);
    refParetoMatched = isParetoRef(idxRef);
    zOrigMatched = double(T_orig.z_eval(idxOrig));
    promotedMask = ~origParetoMatched & refParetoMatched;
    promotedZlt1Mask = promotedMask & isfinite(zOrigMatched) & (zOrigMatched < 1 - 1e-12);
    promotedIdxOrig = idxOrig(promotedZlt1Mask);
    promotedIdxRef = idxRef(promotedZlt1Mask);
    promotedRunKey = matchedRunKey(promotedZlt1Mask);
    promotedTs = matchedTs(promotedZlt1Mask);
    matchedOrigParetoCount = nnz(origParetoMatched);
    matchedRefParetoCount = nnz(refParetoMatched);

    fig = figure("Color", "w", "Toolbar", "none", "Name", "Pareto Frontier Change");
    tiledlayout(fig, 1, 2, "Padding", "compact", "TileSpacing", "compact");
    set(fig, "Position", [80 80 1400 560]);

    % ---- Left: combined frontier (RESULTSSANDBOX figure-5 style) ----
    axL = nexttile; hold(axL, "on");
    scatter(axL, double(Tall.SSdU), double(Tall.SSE), 18, ...
        "filled", "MarkerFaceColor", "k", "MarkerEdgeColor", "none");
    plot_pareto_polyline_with_markers(axL, double(Tf.SSdU), double(Tf.SSE), plotColors(3, :), 2.0, 8, 2.0);
    scatter(axL, double(T1.SSdU(isLeft1)), double(T1.SSE(isLeft1)), 80, plotColors(1,:), ...
        "o", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:), "LineWidth", 1.4);
    scatter(axL, double(T2.SSdU(isLeft2)), double(T2.SSE(isLeft2)), 90, plotColors(2,:), ...
        "^", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:), "LineWidth", 1.4);

    % ---- Right: refinement change ----
    axR = nexttile; hold(axR, "on");
    Torig1 = T_orig(T_orig.run_key == "run1", :);
    Torig2 = T_orig(T_orig.run_key == "run2", :);
    scatter(axR, double(Torig1.SSdU), double(Torig1.SSE), 16, ...
        "o", "filled", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:));
    scatter(axR, double(Torig2.SSdU), double(Torig2.SSE), 20, ...
        "^", "filled", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:));
    scatter(axR, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
        "filled", "MarkerFaceColor", colRefAll, "MarkerEdgeColor", "none");
    TrefPareto = T_ref(isParetoRef, :);
    scatter(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), 80, ...
        "d", "MarkerFaceColor", NATURE_COLOR.ReddishPurple, ...
        "MarkerEdgeColor", NATURE_COLOR.ReddishPurple, "LineWidth", 1.0);

    if ~isempty(promotedIdxRef)
        scatter(axR, double(T_ref.SSdU(promotedIdxRef)), double(T_ref.SSE(promotedIdxRef)), 140, ...
            "p", "MarkerFaceColor", colPromoted, "MarkerEdgeColor", "k", "LineWidth", 0.7);
    end

    for ax = [axL, axR]
        set(ax, "XScale", "log", "YScale", "log", "FontSize", 16);
        xlim(ax, [1e-2, 2e0]);
        ylim(ax, [1e4, 3.5e4]);
        xlabel(ax, "$J_{\mathrm{TV}}$", "Interpreter", "latex");
        ylabel(ax, "$J_{\mathrm{track}}$", "Interpreter", "latex");
        grid(ax, "off");
        box(ax, "off");
        axes(ax);
        format_tick(1, 1);
    end
    title(axL, "$\mathbf{a}$", "Interpreter", "latex");
    title(axR, "$\mathbf{b}$", "Interpreter", "latex");
    axL.TitleHorizontalAlignment = "left";
    axR.TitleHorizontalAlignment = "left";

    % (panel b) no on-figure text annotation by request.

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
    write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef);

    fprintf("Saved: %s\n", outStem + ".png");
    fprintf("Saved: %s\n", outStem + ".pdf");
    fprintf("Saved: %s\n", reportPath);
    fprintf("J_TV rows compared: %d\n", height(T_jtv));
    fprintf("Matched points (run_key + timestamp): %d\n", numel(idxOrig));
    fprintf("Left-panel matched-frontier pool size: %d\n", height(TleftFrontPool));
    fprintf("Left-panel frontier points (non-DOE): %d\n", height(Tf));
    fprintf("Original Pareto points (all original rows after DOE filter): %d\n", nnz(isParetoOrig));
    fprintf("Refined Pareto points (all refined rows): %d\n", nnz(isParetoRef));
    fprintf("Original Pareto points in matched set: %d\n", matchedOrigParetoCount);
    fprintf("Refined Pareto points in matched set: %d\n", matchedRefParetoCount);
    fprintf("Promoted to Pareto with z < 1: %d\n", numel(promotedIdxRef));
end


function print_runtime_cfg(originalFiles, refinedFiles, doeIterationsPerRun)
%PRINT_RUNTIME_CFG Print reproducibility settings used for this comparison.
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
%COMPUTE_JTV_CHANGE_TABLE Quantify objective drift per matched controller.
    [commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), "stable");
    if isempty(commonKeys)
        error("No matching (run_key, timestamp) points between original and refined tables.");
    end
    [runKey, ts] = split_match_key(commonKeys);
    T = table(runKey, ts, ...
        double(T_orig.SSdU(idxOrig)), ...
        double(T_ref.SSdU(idxRef)), ...
        'VariableNames', {'run_key', 'timestamp', 'JTV_original', 'JTV_refined'});
    T.delta_JTV = T.JTV_refined - T.JTV_original;
    T.delta_pct = 100 * T.delta_JTV ./ T.JTV_original;
    T.abs_delta_JTV = abs(T.delta_JTV);
    T = sortrows(T, "abs_delta_JTV", "descend");
end


function print_top1_cfg(projectRoot, T_jtv)
%PRINT_TOP1_CFG Sanity-check the largest change by printing both controller configs.
    if isempty(T_jtv)
        fprintf("\nNo J_TV rows available for cfg comparison.\n");
        return
    end
    runKey = string(T_jtv.run_key(1));
    ts = string(T_jtv.timestamp(1));

    origMat = find_mat_for_timestamp(projectRoot, runKey, ts, true);
    refMat = find_mat_for_timestamp(projectRoot, runKey, ts, false);
    if strlength(origMat) == 0 || strlength(refMat) == 0
        fprintf("\nTop-1 cfg print skipped (missing MAT files) for %s | %s.\n", runKey, ts);
        return
    end

    Sorig = load_cfg_snapshot(origMat);
    Sref = load_cfg_snapshot(refMat);

    fprintf("\n=== Top-1 point by |delta J_TV|: %s | %s ===\n", runKey, ts);
    fprintf("\n--- Original cfg (%s) ---\n", origMat);
    print_cfg_struct(Sorig);

    fprintf("\n--- Refined cfg (%s) ---\n", refMat);
    print_cfg_struct(Sref);
end


function matPath = find_mat_for_timestamp(projectRoot, runKey, ts, isOriginal)
%FIND_MAT_FOR_TIMESTAMP Resolve source artifact for one matched controller identity.
    if isOriginal
        cands = fullfile(projectRoot, "results", runKey, "out_" + ts + ".mat");
    else
        cands = fullfile(projectRoot, "results", "final_fidelity_same_noise", runKey + "_full_f1_same_noise", "out_full_" + ts + ".mat");
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
%PRINT_CFG_STRUCT Display high-value tuning/noise metadata for auditability.
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
    if isfield(S, "cfg_run")
        cfgRun = S.cfg_run;
        if isfield(cfgRun, "mode"), fprintf("cfg_run.mode = %s\n", string(cfgRun.mode)); end
        if isfield(cfgRun, "source_root"), fprintf("cfg_run.source_root = %s\n", string(cfgRun.source_root)); end
        if isfield(cfgRun, "output_root"), fprintf("cfg_run.output_root = %s\n", string(cfgRun.output_root)); end
        if isfield(cfgRun, "results_csv"), fprintf("cfg_run.results_csv = %s\n", string(cfgRun.results_csv)); end
        if isfield(cfgRun, "out_dir"), fprintf("cfg_run.out_dir = %s\n", string(cfgRun.out_dir)); end
        if isfield(cfgRun, "NumWorkers"), fprintf("cfg_run.NumWorkers = %d\n", double(cfgRun.NumWorkers)); end
    end
end


function T = load_and_stack_csv(fileDefs, dropFirstNRows)
%LOAD_AND_STACK_CSV Standardize raw CSV tables into one comparable schema.
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
        if ~isfield(fileDefs(k), "runKey")
            error("fileDefs entry %d is missing required field 'runKey'.", k);
        end
        Tk.run_label = repmat(string(fileDefs(k).runLabel), height(Tk), 1);
        Tk.run_key = repmat(string(fileDefs(k).runKey), height(Tk), 1);
        Tk.iteration = (1:height(Tk)).';
        if dropFirstNRows > 0
            Tk = Tk(Tk.iteration > dropFirstNRows, :);
        end
        if ismember("theta_1", string(Tk.Properties.VariableNames))
            Tk.z_eval = double(Tk.theta_1);
        else
            Tk.z_eval = nan(height(Tk), 1);
        end
        Tk.match_key = Tk.run_key + "|" + string(Tk.timestamp);
        Tk = Tk(:, ["run_label", "run_key", "match_key", "iteration", "timestamp", "SSE", "SSdU", "z_eval"]);
        T = [T; Tk]; %#ok<AGROW>
    end
end


function h = plot_pareto_polyline_with_markers(ax, x, y, colorVal, lineWidth, markerSize, markerLineWidth)
%PLOT_PARETO_POLYLINE_WITH_MARKERS Draw frontier as one visual object (line+markers).
    if isempty(x) || isempty(y)
        h = gobjects(0);
        return
    end
    [xSort, ord] = sort(x(:), "ascend");
    ySort = y(ord);
    h = plot(ax, xSort, ySort, "-o", ...
        "Color", colorVal, ...
        "LineWidth", lineWidth, ...
        "MarkerSize", 15, ...
        "MarkerFaceColor", "w", ...
        "MarkerEdgeColor", colorVal);
    %#ok<INUSD> markerLineWidth retained for call-site compatibility.
end


function h = plot_pareto_continuum(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM Draw a smooth monotone guide through frontier points.
    [xSort, ord] = sort(x(:), "ascend");
    ySort = y(ord);

    if numel(xSort) < 3
        h = plot(ax, xSort, ySort, "-", "Color", curveColor, "LineWidth", 2.0);
        return;
    end

    lx = log10(xSort);
    ly = log10(ySort);
    lxqIn = linspace(min(lx), max(lx), 220);
    lyqIn = pchip(lx, ly, lxqIn);
    xq = (10.^lxqIn).';
    yq = (10.^lyqIn).';

    if nargin >= 6 && ~isempty(yBounds)
        yTop = max(yBounds);
        if yTop > yq(1)
            xq = [xq(1); xq];
            yq = [yTop; yq];
        end
    end

    if nargin >= 5 && ~isempty(xBounds)
        xRight = max(xBounds);
        if xRight > xq(end)
            xq = [xq; xRight];
            yq = [yq; yq(end)];
        end
    end

    h = plot(ax, xq, yq, "-", "Color", curveColor, "LineWidth", 2.0);
end


function [TrefKeep, info] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig, T_orig_full)
%KEEP_REFINED_FROM_ORIGINAL_PARETO Enforce refined eligibility for fair comparison.
    refKeys = string(T_ref.match_key);
    origKeys = string(T_orig.match_key);
    origFullKeys = string(T_orig_full.match_key);
    origParetoKeys = string(T_orig.match_key(isParetoOrig));

    inOrigFull = ismember(refKeys, origFullKeys);
    inOrig = ismember(refKeys, origKeys);
    inOrigPareto = ismember(refKeys, origParetoKeys);

    info = struct();
    info.n_ref_input = height(T_ref);
    info.n_keep_before_dedup = nnz(inOrigPareto);
    info.n_drop_doe = nnz(inOrigFull & ~inOrig);
    info.n_drop_not_in_orig = nnz(~inOrigFull);
    info.n_drop_in_orig_not_pareto = nnz(inOrig & ~inOrigPareto);
    info.n_drop_total = info.n_drop_doe + info.n_drop_not_in_orig + info.n_drop_in_orig_not_pareto;
    info.dropped_doe_keys = refKeys(inOrigFull & ~inOrig);

    TrefKeep = T_ref(inOrigPareto, :);

    % Guard against repeated refined evaluations of the same key.
    if ~isempty(TrefKeep)
        [~, uniqIdx] = unique(string(TrefKeep.match_key), "stable");
        info.n_drop_duplicate_refined = height(TrefKeep) - numel(uniqIdx);
        TrefKeep = TrefKeep(uniqIdx, :);
    else
        info.n_drop_duplicate_refined = 0;
    end
    info.n_keep_final = height(TrefKeep);
end


function print_refined_filter_summary(info)
%PRINT_REFINED_FILTER_SUMMARY Report why refined points were kept or dropped.
    fprintf("\n=== Refined sample eligibility check ===\n");
    fprintf("Input refined rows: %d\n", info.n_ref_input);
    fprintf("Dropped (mapped to original DOE rows): %d\n", info.n_drop_doe);
    fprintf("Dropped (not found in original run datasets): %d\n", info.n_drop_not_in_orig);
    fprintf("Dropped (original but not Pareto): %d\n", info.n_drop_in_orig_not_pareto);
    fprintf("Dropped duplicate refined rows: %d\n", info.n_drop_duplicate_refined);
    fprintf("Kept refined rows (original non-DOE Pareto-linked): %d\n", info.n_keep_final);
    if info.n_drop_doe > 0
        fprintf("Dropped DOE-linked keys:\n");
        for i = 1:numel(info.dropped_doe_keys)
            fprintf("  - %s\n", string(info.dropped_doe_keys(i)));
        end
    end
    if info.n_drop_total == 0 && info.n_drop_duplicate_refined == 0
        fprintf("Verified: all refined samples map to original non-DOE Pareto points.\n");
    end
end


function draw_frontier(ax, Tpareto, colorVal, legendLabel)
%DRAW_FRONTIER Legacy helper for frontier overlays (kept for compatibility).
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
%COMPUTE_PARETO_MASK Identify non-dominated tradeoffs (min SSE, min SSdU).
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


function write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef)
%WRITE_PROMOTED_REPORT Persist promoted-point audit rows for manuscript traceability.
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
    fprintf(fid, "run_key,timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n");

    for i = 1:numel(promotedTs)
        io = promotedIdxOrig(i);
        ir = promotedIdxRef(i);
        fprintf(fid, "%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n", ...
            promotedRunKey(i), promotedTs(i), ...
            double(T_orig.SSE(io)), double(T_orig.SSdU(io)), double(T_orig.z_eval(io)), ...
            double(T_ref.SSE(ir)), double(T_ref.SSdU(ir)));
    end
end


function S = load_cfg_snapshot(matPath)
%LOAD_CFG_SNAPSHOT Load only audit-relevant MAT fields to avoid handle issues.
    S = struct();
    vars = string(who("-file", matPath));
    if any(vars == "out")
        tmp = load(matPath, "out");
        S.out = tmp.out;
    end
    if any(vars == "cfg_run")
        tmp = load(matPath, "cfg_run");
        S.cfg_run = tmp.cfg_run;
    end
end


function [runKey, ts] = split_match_key(matchKey)
%SPLIT_MATCH_KEY Decode controller identity keys back to run and timestamp.
    n = numel(matchKey);
    runKey = strings(n, 1);
    ts = strings(n, 1);
    for i = 1:n
        parts = split(string(matchKey(i)), "|");
        if numel(parts) >= 2
            runKey(i) = parts(1);
            ts(i) = parts(2);
        else
            runKey(i) = "";
            ts(i) = string(matchKey(i));
        end
    end
end
