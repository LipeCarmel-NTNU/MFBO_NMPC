function timestamps = load_results_csv_timestamps(results_csv)
%LOAD_RESULTS_CSV_TIMESTAMPS Read the timestamp column of a results CSV.
%
%   Returns an empty string array when the file is absent or has no
%   timestamp column.

    if ~isfile(results_csv)
        timestamps = strings(0, 1);
        return
    end
    T = readtable(results_csv, "TextType", "string");
    if ~ismember("timestamp", string(T.Properties.VariableNames))
        timestamps = strings(0, 1);
        return
    end
    timestamps = unique(string(T.timestamp), "stable");
end
