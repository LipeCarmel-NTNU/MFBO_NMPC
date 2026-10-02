function test_realtime()
%TEST_REALTIME Edits 1 to 4, MATLAB side: the real-time feasibility deadline.

    repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repo));
    configure_pool(false, 0);
    warning("off", "NMPC:fallback");
    warning("off", "NMPC:failed_again");

    rng(123);
    base = nmpc_base(sigma_y = [0.001 0.1 0.1]);
    theta = [0.02, 0, 1, 0 0 0, -1000 -1000 -1000, 0 0 0];
    quiet = {"horizon", "fidelity", "extrapolate", false, ...
             "terminal_cost", "lqr", "verbosity", "quiet"};

    %% The value: 60 s step times the 1.10 tolerance
    assert(abs(base.rt_deadline_s - 66) < 1e-9, "deadline is %g s, not 66", base.rt_deadline_s);
    assert(base.rt_tolerance == 1.10);
    fprintf("deadline value: PASS (%g s)\n", base.rt_deadline_s);

    %% The controller carries it, and Inf is the off state of a bare NMPC
    cfg = decode_theta(theta, base.nx, base.nu);
    nmpc = build_nmpc(base, cfg, terminal_cost = "lqr");
    assert(nmpc.rt_deadline_s == base.rt_deadline_s, "build_nmpc did not copy the deadline");
    probe = NMPC_defaults();
    assert(isinf(probe), "NMPC default rt_deadline_s is %g, not Inf", probe);
    fprintf("build_nmpc wiring: PASS\n");

    %% The output function stops fmincon once the budget is spent
    nmpc.rt_deadline_s = 1e-6;
    nmpc.rt_t0 = tic;
    nmpc.rt_violated = false;
    pause(0.01);
    x0 = [1.0 10 0];
    nmpc.solve(x0, zeros(1, base.nu));
    assert(nmpc.rt_violated, "rt_violated not set by the output function");
    assert(nmpc.latest_flag == -1, "fmincon exit flag %g, expected -1 from the stop", nmpc.latest_flag);
    assert(isempty(nmpc.latest_wopt), "a stopped solve poisoned the warm start");
    fprintf("output function stop: PASS\n");

    %% Disarmed, the same controller solves normally
    nmpc.rt_deadline_s = Inf;
    nmpc.rt_violated = false;
    nmpc.solve(x0, zeros(1, base.nu));
    assert(~nmpc.rt_violated && nmpc.latest_flag >= 0, "a disarmed solve was stopped");
    fprintf("disarmed solve: PASS\n");

    %% A whole evaluation under the real 66 s deadline completes
    out = simulate_nmpc(base, theta, quiet{:}, run_id = "rt_ok");
    assert(isfinite(out.SSE));
    fprintf("66 s deadline, cheap theta: PASS (SSE=%.6g)\n", out.SSE);

    %% A deadline no step can meet aborts the evaluation with the identifier
    tight = base;
    tight.rt_deadline_s = 1e-3;
    ckpt = string(tempname) + ".mat";
    cleanup = onCleanup(@() delete_if_exists(ckpt));
    try
        simulate_nmpc(tight, theta, quiet{:}, run_id = "rt_bad", ...
            checkpoint_path = ckpt, checkpoint_id = "eval_1");
        error("test:noAbort", "the evaluation completed under a 1 ms deadline");
    catch ME
        assert(ME.identifier == "NMPC:realtimeInfeasible", "wrong error: %s (%s)", ME.identifier, ME.message);
        fprintf("tight deadline aborts: PASS (%s)\n", ME.message);
    end
    assert(~isfile(ckpt), "a checkpoint was written for the aborted step");
    fprintf("no checkpoint for the aborted step: PASS\n");

    %% serve_requests does not count the abort toward the failure limit
    src = fileread(fullfile(repo, "dependencies", "io", "serve_requests.m"));
    assert(contains(src, 'ME.identifier == "NMPC:realtimeInfeasible"'), ...
        "serve_requests does not single out the real-time abort");
    fprintf("serve_requests guard: PASS\n");
end

function v = NMPC_defaults()
    mc = ?NMPC;
    p = mc.PropertyList(strcmp({mc.PropertyList.Name}, "rt_deadline_s"));
    assert(~isempty(p), "NMPC has no rt_deadline_s property");
    v = p.DefaultValue;
end
