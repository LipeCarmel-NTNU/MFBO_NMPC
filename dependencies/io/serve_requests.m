function n_served = serve_requests(cfg, evaluate_fn)
%SERVE_REQUESTS Poll the inbox and evaluate one request at a time.
%
%   The driver process writes an evaluation request to cfg.theta_txt and waits
%   for the matching row to appear in the results CSV. This function is the
%   MATLAB half of that exchange. It is shared by the initialisation and the
%   optimisation entry points so that both handle interruption the same way.
%
%   evaluate_fn is called as evaluate_fn(req) with the request struct returned
%   by read_theta_from_txt. It is expected to append exactly one results row
%   carrying req.eval_id.
%
%   An error raised by evaluate_fn does not end the loop. The failure is
%   appended to cfg.failures_csv, which the driver polls alongside the results
%   file, and serving continues. A run of consecutive failures larger than
%   cfg.max_consecutive_failures does end the loop, because past that point the
%   fault is in the configuration rather than in one candidate and continuing
%   would burn the budget on rows that never reach the optimiser.
%
%   Required cfg fields:
%     theta_txt, results_csv, failures_csv, log_path, theta_len, poll_s,
%     lock_path, lock_stale_s, max_consecutive_failures
%   Optional:
%     stop_after   return once this many requests have been served

    if ~isfield(cfg, "stop_after")
        cfg.stop_after = Inf;
    end

    n_served = 0;
    n_consecutive_failures = 0;
    last_signature = "";
    last_eval_id = -Inf;

    fprintf("Serving requests from %s (Ctrl-C to stop).\n", cfg.theta_txt);

    while n_served < cfg.stop_after
        [req, ok] = read_theta_from_txt(cfg.theta_txt, cfg.theta_len);

        if ~ok
            pause(cfg.poll_s);
            continue
        end

        if req.signature == last_signature
            pause(cfg.poll_s);
            continue
        end

        % A request whose index is not ahead of the last one served is a file
        % this process has already answered, which is what a restarted driver
        % leaves behind. Serving it again would append a duplicate row under an
        % eval_id the driver has already recorded.
        if req.eval_id <= last_eval_id
            last_signature = req.signature;
            pause(cfg.poll_s);
            continue
        end

        acquire_lock(cfg.lock_path, cfg.lock_stale_s);
        cleanup_lock = onCleanup(@() delete_if_exists(cfg.lock_path)); %#ok<NASGU>

        t_start = tic;
        try
            evaluate_fn(req);
            n_consecutive_failures = 0;
            fprintf("  eval %d served in %.1f s\n", req.eval_id, toc(t_start));
        catch ME
            n_consecutive_failures = n_consecutive_failures + 1;
            append_failure_row(cfg.failures_csv, req.eval_id, ...
                timestamp_compact(), ME.identifier, ME.message, req.theta);
            if strlength(cfg.log_path) > 0
                log_simulation_event(cfg.log_path, sprintf( ...
                    "eval %d failed (%d in a row): %s | %s", ...
                    req.eval_id, n_consecutive_failures, ME.identifier, ME.message));
            end
            fprintf(2, "  eval %d FAILED (%d in a row): %s\n", ...
                req.eval_id, n_consecutive_failures, ME.message);
            fprintf(2, "%s\n", getReport(ME, "extended", "hyperlinks", "off"));
        end

        clear cleanup_lock
        delete_if_exists(cfg.lock_path);

        last_signature = req.signature;
        last_eval_id = req.eval_id;
        n_served = n_served + 1;

        if n_consecutive_failures > cfg.max_consecutive_failures
            error("serve_requests:tooManyFailures", ...
                ["%d consecutive evaluations failed, which is past the limit of " ...
                 "%d. The fault is unlikely to be candidate-specific. See %s."], ...
                n_consecutive_failures, cfg.max_consecutive_failures, cfg.failures_csv);
        end
    end
end


function acquire_lock(lock_path, stale_s)
%ACQUIRE_LOCK Create the busy marker, clearing one left behind by a crash.
% The driver waits on this file before writing a new request. A process that
% died mid-evaluation leaves it in place, so a lock older than stale_s is
% treated as abandoned and removed with a warning rather than trusted.
    if isfile(lock_path)
        d = dir(lock_path);
        age_s = seconds(datetime("now") - datetime(d.datenum, "ConvertFrom", "datenum"));
        if age_s > stale_s
            warning("serve_requests:staleLock", ...
                "Removing a lock file %.0f s old, left by an earlier process: %s", ...
                age_s, lock_path);
        end
        delete_if_exists(lock_path);
    end

    fid = fopen(lock_path, "w");
    if fid < 0
        warning("serve_requests:lockFailed", ...
            "Could not create the lock file (pwd = %s): %s", pwd, lock_path);
        return
    end
    fprintf(fid, "pid unknown, created %s\n", timestamp_compact());
    fclose(fid);
end
