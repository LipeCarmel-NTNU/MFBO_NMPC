function cfg_run = server_config(serves_phase)
%SERVER_CONFIG Settings that both MATLAB servers share.
%
%   cfg_run = server_config(serves_phase) returns the part of cfg_run that
%   main_initialization (serves_phase 0) and main_BO (serves_phase 1) have in
%   common, and results_root, the folder of the campaign. The caller adds the
%   paths that differ between the phases and theta_len, which depends on base.
%
%   The results folder comes from MFBO_RESULTS_DIR. The driver derives it from
%   the case and the Sobol seed (RESULTS_DIR in pipeline/matlab_interface.py),
%   and run_supervised.py passes it in that variable. A server does not know the
%   case, so it cannot derive the folder itself. The exchange stays at the
%   project root.

    arguments
        serves_phase (1,1) double {mustBeMember(serves_phase, [0 1])}
    end

    results_root = string(getenv("MFBO_RESULTS_DIR"));
    if strlength(results_root) == 0
        error("MFBO:resultsDir", ...
            "MFBO_RESULTS_DIR is not set. run_supervised.py sets it. When you " + ...
            "start MATLAB yourself, first run setenv(""MFBO_RESULTS_DIR"", " + ...
            """<folder>"") with the folder that run_pipeline.py prints.");
    end

    cfg_run = struct();
    cfg_run.results_root = results_root;
    cfg_run.theta_txt = fullfile("inbox", "theta.txt");
    cfg_run.poll_s = 2.0;
    % One failures file for both phases, because the serve loop does not know
    % the phase of a request that raised.
    cfg_run.failures_csv = fullfile(results_root, "failures.csv");
    cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
    cfg_run.lock_path = "matlab.lock";
    cfg_run.lock_stale_s = 6 * 3600;
    % The live value. max_consecutive_failures in run_config.py only records it
    % in the manifest.
    cfg_run.max_consecutive_failures = 5;
    % run_idun.slurm reserves one core more than this, in MATLAB_WORKERS.
    cfg_run.NumWorkers = 8;
    cfg_run.serves_phase = serves_phase;
end
