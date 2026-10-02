function test_server_config()
%TEST_SERVER_CONFIG The settings both servers share live in one function.
%   The expected values are the literals main_initialization.m and main_BO.m
%   set before the extraction.

    repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repo));
    root = "results/MF_case2_s123";
    saved = getenv("MFBO_RESULTS_DIR");
    restore = onCleanup(@() setenv("MFBO_RESULTS_DIR", saved)); %#ok<NASGU>
    setenv("MFBO_RESULTS_DIR", root);

    for phase = [0 1]
        c = server_config(phase);
        assert(c.results_root == root);
        assert(isequal(c.theta_txt, fullfile("inbox", "theta.txt")));
        assert(c.poll_s == 2.0);
        assert(isequal(c.failures_csv, fullfile(root, "failures.csv")));
        assert(isequal(c.log_path, fullfile("SIMULATIONS_LOG.txt")));
        assert(c.lock_path == "matlab.lock");
        assert(c.lock_stale_s == 10 * 3600);
        assert(c.max_consecutive_failures == 5);
        assert(c.NumWorkers == 8);
        assert(c.serves_phase == phase);
        fprintf("server_config(%d): PASS\n", phase);
    end

    setenv("MFBO_RESULTS_DIR", "");
    try
        server_config(0);
        error("test:noError", "no error with MFBO_RESULTS_DIR unset");
    catch ME
        assert(ME.identifier == "MFBO:resultsDir", "wrong error %s", ME.identifier);
    end
    fprintf("unset variable: PASS\n");

    % The servers keep only what differs between the phases.
    shared = ["poll_s", "lock_stale_s", "lock_path", "max_consecutive_failures", ...
              "theta_txt", "failures_csv", "log_path", "sigma_y", "NumWorkers"];
    for m = ["main_initialization.m", "main_BO.m"]
        src = fileread(fullfile(repo, m));
        set_here = @(f) ~isempty(regexp(src, "cfg_run\." + f + "\s*=", "once"));
        left = shared(arrayfun(set_here, shared));
        for extra = ["getenv(""MFBO_RESULTS_DIR"")", "NumWorkers = ", "rng("]
            if contains(src, extra), left(end+1) = extra; end %#ok<AGROW>
        end
        assert(isempty(left), "%s still sets: %s", m, strjoin(left, ", "));
        assert(contains(src, "server_config(" + double(m == "main_BO.m") + ")"), "%s does not call server_config", m);
        fprintf("%s: PASS\n", m);
    end
end
