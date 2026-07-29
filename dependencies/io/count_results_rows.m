function [n, last_eval_id] = count_results_rows(results_csv)
%COUNT_RESULTS_ROWS Number of data rows in a results CSV, and the last eval_id.
%
%   Counts lines rather than parsing the file, because the caller only needs to
%   know how much of the budget is already spent when resuming. A trailing line
%   without a newline still counts, and the header does not.
%
%   last_eval_id is the first field of the final data row, or 0 when the file
%   holds none. Resuming from it rather than from the row count is what keeps
%   the indices correct after a failed evaluation, which produces no results
%   row but does consume an eval_id.

    n = 0;
    last_eval_id = 0;

    if ~isfile(results_csv)
        return
    end

    fid = fopen(results_csv, "r");
    if fid < 0
        error("Could not open results file for reading: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    header = fgetl(fid);
    if ~ischar(header)
        return
    end

    last_line = "";
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if strlength(strtrim(string(line))) == 0
            continue
        end
        n = n + 1;
        last_line = string(line);
    end

    if n > 0
        parts = split(last_line, ",");
        v = str2double(parts(1));
        if isfinite(v)
            last_eval_id = v;
        end
    end
end
