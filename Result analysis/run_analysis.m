function run_analysis()
%RUN_ANALYSIS Execute the full analysis pipeline and sync PDFs to Overleaf.
%
% Pipeline (in order):
% 1) Result analysis/analyze_test_run_metrics.m
% 2) Result analysis/recover_boxplot_data.py
% 3) Result analysis/resultssandbox.m
% 4) J surrogate/surrogate.m
% 5) results/numerical results/compile_results_txt.m
% 6) Copy and overwrite all PDFs from results/graphical_results to
%    Overleaf---MFBO-NMPC.

    scriptDir = fileparts(mfilename("fullpath"));
    projectRoot = fileparts(scriptDir);
    addpath(genpath(fullfile(projectRoot, "dependencies")));

    fprintf("\n=== RUN_ANALYSIS: starting full pipeline ===\n");

    run_script_in_base(fullfile(scriptDir, "analyze_test_run_metrics.m"), ...
        "analyze_test_run_metrics");

    run_python_recover_boxplot(fullfile(scriptDir, "recover_boxplot_data.py"));

    run_script_in_base(fullfile(scriptDir, "resultssandbox.m"), ...
        "resultssandbox");

    run_script_in_base(fullfile(projectRoot, "J surrogate", "surrogate.m"), ...
        "surrogate");

    run_script_in_base(fullfile(projectRoot, "results", "numerical results", "compile_results_txt.m"), ...
        "compile_results_txt");

    copy_graphical_pdfs_to_overleaf( ...
        fullfile(projectRoot, "results", "graphical_results"), ...
        fullfile(projectRoot, "Overleaf---MFBO-NMPC"));

    fprintf("=== RUN_ANALYSIS: completed successfully ===\n\n");
end


function run_script_in_base(scriptPath, label)
%RUN_SCRIPT_IN_BASE Execute a MATLAB script path in base workspace.
    if ~isfile(scriptPath)
        error("Script not found: %s", scriptPath);
    end
    fprintf("\n[RUN] %s\n", label);
    cmd = sprintf("run('%s');", strrep(scriptPath, "'", "''"));
    evalin("base", cmd);
end


function run_python_recover_boxplot(pyPath)
%RUN_PYTHON_RECOVER_BOXPLOT Run recover_boxplot_data.py via system Python.
    if ~isfile(pyPath)
        error("Python script not found: %s", pyPath);
    end

    fprintf("\n[RUN] recover_boxplot_data.py\n");
    pyExec = detect_python_executable();
    if strlength(pyExec) == 0
        warning("No Python executable found. Skipping recover_boxplot_data.py.");
        return
    end

    cmd = sprintf('%s "%s"', pyExec, pyPath);
    status = system(cmd);
    if status ~= 0
        warning("Python script failed with status %d: %s", status, pyPath);
    end
end


function pyExec = detect_python_executable()
%DETECT_PYTHON_EXECUTABLE Return a working Python command string.
    if ispc
        candidates = ["python", "py -3", "py"];
    else
        candidates = ["python3", "python"];
    end

    pyExec = "";
    for i = 1:numel(candidates)
        testCmd = sprintf("%s --version", candidates(i));
        status = system(testCmd);
        if status == 0
            pyExec = candidates(i);
            return
        end
    end
end


function copy_graphical_pdfs_to_overleaf(srcDir, dstDir)
%COPY_GRAPHICAL_PDFS_TO_OVERLEAF Copy and overwrite generated PDFs to Overleaf.
    fprintf("\n[SYNC] Copy PDFs to Overleaf (overwrite enabled)\n");

    if ~isfolder(srcDir)
        error("Source folder not found: %s", srcDir);
    end
    if ~isfolder(dstDir)
        error("Overleaf folder not found: %s", dstDir);
    end

    pdfFiles = dir(fullfile(srcDir, "*.pdf"));
    if isempty(pdfFiles)
        warning("No PDFs found in %s", srcDir);
        return
    end

    for i = 1:numel(pdfFiles)
        srcPath = fullfile(pdfFiles(i).folder, pdfFiles(i).name);
        dstPath = fullfile(dstDir, pdfFiles(i).name);
        ok = copyfile(srcPath, dstPath, "f"); % force overwrite
        if ~ok
            error("Failed to copy %s to %s", srcPath, dstPath);
        end
        fprintf("  copied: %s\n", pdfFiles(i).name);
    end
end
