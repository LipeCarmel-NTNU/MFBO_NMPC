function init_results_csv(results_csv, theta_len)
%INIT_RESULTS_CSV Write the results CSV header if the file does not exist.
%
%   Columns are timestamp, SSE, SSdU, J, runtime_s followed by theta_1 to
%   theta_<theta_len>. An existing file is left untouched so that a run can
%   append to it after an interruption.

    if isfile(results_csv)
        return
    end

    fid = fopen(results_csv, "w");
    if fid < 0
        error("Could not open results file for writing: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "timestamp,SSE,SSdU,J,runtime_s");
    for k = 1:theta_len
        fprintf(fid, ",theta_%d", k);
    end
    fprintf(fid, "\n");
end
