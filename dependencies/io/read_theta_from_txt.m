function [theta, signature, ok] = read_theta_from_txt(theta_txt)
%READ_THETA_FROM_TXT Read theta from the last non-empty line of a text file.
%
%   ok is false when the file is missing, unreadable, empty or unparseable, so
%   the caller can poll without a try/catch.
%
%   signature combines the file modification time with the content of that
%   line. It changes only when a new theta arrives, which lets the caller skip
%   a file it has already evaluated.

    theta = [];
    signature = "";
    ok = false;

    if ~isfile(theta_txt)
        return
    end

    try
        raw = fileread(theta_txt);
    catch
        return
    end

    lines = splitlines(string(raw));
    lines = strip(lines);
    lines = lines(lines ~= "");

    if isempty(lines)
        return
    end

    last_line = lines(end);
    vals = sscanf(last_line, "%f").';
    if isempty(vals)
        return
    end

    theta = vals;

    d = dir(theta_txt);
    if isempty(d)
        signature = "missing";
    else
        signature = string(d.datenum) + "|" + string(last_line);
    end
    ok = true;
end
