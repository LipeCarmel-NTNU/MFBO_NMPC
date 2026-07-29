function check_results_header(results_csv, theta_len)
%CHECK_RESULTS_HEADER Refuse to append to a CSV written under another schema.
%
%   A run that resumes into a file whose columns differ from the current
%   definition would append rows that parse without error and mean something
%   else. The header is compared verbatim against results_csv_header and any
%   difference stops the run, naming both versions.
%
%   A missing file is not an error: init_results_csv creates it.

    if ~isfile(results_csv)
        return
    end

    fid = fopen(results_csv, "r");
    if fid < 0
        error("Could not open results file for reading: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    first = fgetl(fid);
    if ~ischar(first)
        error("Results file exists but is empty: %s", results_csv);
    end

    expected = results_csv_header(theta_len);
    if strtrim(string(first)) ~= expected
        error(['Results file was written under a different schema: %s\n' ...
               '  found:    %s\n' ...
               '  expected: %s\n' ...
               'Move the old file aside before resuming, or point the run at a ' ...
               'new results directory.'], results_csv, strtrim(string(first)), expected);
    end
end
