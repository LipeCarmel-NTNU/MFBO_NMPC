function rerun_stored_controllers(cfg_run, base, run_label, timestamps, opts)
%RERUN_STORED_CONTROLLERS Rerun stored controllers at full horizon.
%
%   Loads out_<timestamp>.mat from <cfg_run.source_root>/<run_label>, forces
%   theta(1) = 1 so the cost is measured over the whole horizon, simulates, and
%   writes one summary row plus one output file per controller into
%   <cfg_run.output_root>/<run_label><out_suffix>.
%
%   A controller whose output file already exists is skipped, so an interrupted
%   sweep resumes where it stopped.
%
%   Name-value options:
%     out_suffix       appended to run_label to form the output folder name
%     recover_csv_rows when a result exists on disk but has no CSV row, rebuild
%                      the row from the stored struct instead of resimulating
%     write_list       record the selected timestamps next to the results

    arguments
        cfg_run struct
        base struct
        run_label (1,1) string
        timestamps (:,1) string
        opts.out_suffix (1,1) string = "_full"
        opts.recover_csv_rows (1,1) logical = false
        opts.write_list (1,1) logical = false
    end

    src_dir = fullfile(cfg_run.source_root, run_label);
    out_dir = fullfile(cfg_run.output_root, run_label + opts.out_suffix);
    ensure_dir(out_dir);

    if opts.write_list
        write_timestamp_list(out_dir, timestamps);
    end

    results_csv = fullfile(out_dir, "results_full.csv");
    theta_len = 1 + 2 + base.nx + 2*base.nu;
    check_results_header(results_csv, theta_len);
    init_results_csv(results_csv, theta_len);
    existing_ts = load_results_csv_timestamps(results_csv);

    for i = 1:numel(timestamps)
        ts = timestamps(i);
        mat_path = fullfile(src_dir, "out_" + ts + ".mat");
        if ~isfile(mat_path)
            warning("Missing controller file: %s", mat_path);
            continue
        end

        full_result_path = fullfile(out_dir, "out_full_" + ts + ".mat");
        if isfile(full_result_path)
            if opts.recover_csv_rows && ~ismember(ts, existing_ts)
                if recover_csv_row(results_csv, full_result_path, ts, i)
                    existing_ts(end+1,1) = ts; %#ok<AGROW>
                end
            else
                fprintf("Skipping %s (results already exist in %s)\n", ts, full_result_path);
            end
            continue
        end

        S = load(mat_path, "out");
        if ~isfield(S, "out")
            warning("No 'out' struct in %s", mat_path);
            continue
        end

        theta = S.out.theta;
        theta(1) = 1;                       % full simulation horizon

        out = simulate_nmpc(base, theta, ...
            horizon = "fidelity", ...
            extrapolate = false, ...
            terminal_cost = "lqr", ...
            run_id = ts, ...
            log_path = string(cfg_run.log_path));

        % eval_id is the position in the timestamp list. These reruns are not
        % served through the inbox, so no driver index exists; the position is
        % stable for a given list and keeps the schema uniform.
        append_results_row(results_csv, i, char(ts), "FULL", out, theta);
        existing_ts(end+1,1) = ts; %#ok<AGROW>
        save(full_result_path, "ts", "out", "theta", "cfg_run", "base");
    end
end

function ok = recover_csv_row(results_csv, full_result_path, ts, eval_id)
%RECOVER_CSV_ROW Rebuild a missing summary row from a stored result.
    ok = false;
    R = load(full_result_path, "out", "theta");
    if ~isfield(R, "out") || ~isfield(R, "theta") ...
            || ~isfield(R.out, "SSE") || ~isfield(R.out, "SSdU") || ~isfield(R.out, "runtime_s")
        warning("Existing result is missing required fields; cannot recover CSV row: %s", full_result_path);
        return
    end
    append_results_row(results_csv, eval_id, char(ts), "FULL", R.out, R.theta);
    fprintf("Recovered missing CSV row for %s from existing %s\n", ts, full_result_path);
    ok = true;
end

function write_timestamp_list(out_dir, timestamps)
%WRITE_TIMESTAMP_LIST Record the selected fidelity-1 timestamps.
    path = fullfile(out_dir, "selected_f1_timestamps.txt");
    fid = fopen(path, "w");
    if fid < 0
        warning("Could not open timestamp list for writing: %s", path);
        return
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "Selected timestamps with theta_1 == 1\n");
    fprintf(fid, "Count: %d\n\n", numel(timestamps));
    for i = 1:numel(timestamps)
        fprintf(fid, "%s\n", timestamps(i));
    end
end
