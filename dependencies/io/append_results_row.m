function append_results_row(results_csv, ts, SSE, SSdU, J, runtime_s, theta)
%APPEND_RESULTS_ROW Append one evaluation to the results CSV.
%
%   Values are written with %.17g so that a double survives the round trip
%   through text.

    fid = fopen(results_csv, "a");
    if fid < 0
        error("Could not open results file for appending: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "%s,%.17g,%.17g,%.17g,%.17g", ts, SSE, SSdU, J, runtime_s);

    theta = theta(:).';
    for k = 1:numel(theta)
        fprintf(fid, ",%.17g", theta(k));
    end
    fprintf(fid, "\n");
end
