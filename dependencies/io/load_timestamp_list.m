function timestamps = load_timestamp_list(path)
%LOAD_TIMESTAMP_LIST Read run timestamps from a plain text list.
%
%   Every non-empty line starting with "20" counts as a timestamp, so
%   headers and comment lines in the same file are ignored. Duplicates are
%   removed and the original order is kept.

    if ~isfile(path)
        error("Timestamp file not found: %s", path);
    end
    lines = strip(string(splitlines(fileread(path))));
    is_ts = lines ~= "" & startsWith(lines, "20");
    timestamps = unique(lines(is_ts), "stable");
end
