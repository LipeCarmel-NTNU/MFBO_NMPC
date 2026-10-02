function test_doe_checkpoint()
%TEST_DOE_CHECKPOINT Edit 5: design-phase checkpointing.
%   Static: run_and_log in main_initialization.m passes a checkpoint keyed on
%   eval_id and deletes it after the results row.
%   Behavioural: a DOE-configured evaluation interrupted in case 2 resumes
%   under a fresh run_id and reproduces the uninterrupted cost.

    repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repo));
    configure_pool(false, 0);

    %% Static wiring
    src = fileread(fullfile(repo, "main_initialization.m"));
    must = ["checkpoint_path = ckpt_path", ...
            "checkpoint_id = sprintf(""eval_%d"", req.eval_id)", ...
            "sprintf(""eval_%d.mat"", req.eval_id)", ...
            "delete_if_exists(ckpt_path)"];
    for s = must
        assert(contains(src, s), "main_initialization.m is missing: %s", s);
    end
    i_row = strfind(src, "append_results_row(");
    i_del = strfind(src, "delete_if_exists(ckpt_path)");
    assert(i_del(1) > i_row(1), "checkpoint is deleted before the results row");
    fprintf("static wiring: PASS\n");

    %% Behavioural resume
    rng(123);
    base = nmpc_base(sigma_y = [0.001 0.1 0.1]);
    theta = [0.02, 0, 1, 0 0 0, -1000 -1000 -1000, 0 0 0];
    doe = {"horizon", "fidelity", "extrapolate", false, ...
           "terminal_cost", "lqr", "verbosity", "quiet"};

    ref = simulate_nmpc(base, theta, doe{:}, run_id = "ref");

    ckpt = string(tempname) + ".mat";
    cleanup = onCleanup(@() delete_if_exists(ckpt));
    crash = @(NMPC, t, xk, case_id, sp, i, first) crash_in_case2(sp, case_id, i);
    try
        simulate_nmpc(base, theta, doe{:}, run_id = "attempt1", ...
            checkpoint_path = ckpt, checkpoint_id = "eval_7", setpoint_fn = crash);
        error("test:noCrash", "the injected crash did not fire");
    catch ME
        assert(ME.identifier == "test:crash", "unexpected error: %s", ME.message);
    end
    assert(isfile(ckpt), "no checkpoint left after the crash");

    res = simulate_nmpc(base, theta, doe{:}, run_id = "attempt2", ...
        checkpoint_path = ckpt, checkpoint_id = "eval_7");

    fprintf("ref SSE=%.17g SSdU=%.17g\nres SSE=%.17g SSdU=%.17g\n", ...
        ref.SSE, ref.SSdU, res.SSE, res.SSdU);
    assert(abs(res.SSE - ref.SSE) <= 1e-9 * max(1, abs(ref.SSE)), "SSE differs after resume");
    assert(abs(res.SSdU - ref.SSdU) <= 1e-9 * max(1, abs(ref.SSdU)), "SSdU differs after resume");
    fprintf("behavioural resume: PASS\n");
end

function sp = crash_in_case2(sp, case_id, i)
    if case_id == 2 && i == 3
        error("test:crash", "injected crash");
    end
end
