function append_failure_row(failures_csv, eval_id, ts, identifier, message, theta)
%APPEND_FAILURE_ROW Record an evaluation that raised, so the driver can react.
%
%   The driver polls this file alongside the results CSV. Without it, a failed
%   evaluation is indistinguishable from a slow one and the driver waits for a
%   row that will never arrive.
%
%   The message is flattened to one line and quoted, because an MException
%   report can carry newlines and commas.

    ensure_parent_dir(failures_csv);

    if ~isfile(failures_csv)
        fid = fopen(failures_csv, "w");
        if fid < 0
            error("Could not open failures file for writing: %s", failures_csv);
        end
        c = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid, "eval_id,timestamp,identifier,message");
        for k = 1:numel(theta)
            fprintf(fid, ",theta_%d", k);
        end
        fprintf(fid, "\n");
        clear c
    end

    flat = regexprep(string(message), "\s+", " ");
    flat = replace(flat, """", "''");

    line = sprintf("%d,%s,%s,""%s""", eval_id, ts, string(identifier), flat);
    theta = theta(:).';
    for k = 1:numel(theta)
        line = line + sprintf(",%.17g", theta(k));
    end

    fid = fopen(failures_csv, "a");
    if fid < 0
        error("Could not open failures file for appending: %s", failures_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "%s\n", line);
end
