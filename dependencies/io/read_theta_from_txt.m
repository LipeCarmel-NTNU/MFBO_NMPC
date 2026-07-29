function [req, ok] = read_theta_from_txt(theta_txt, theta_len)
%READ_THETA_FROM_TXT Read one evaluation request from the inbox file.
%
%   The driver writes a single line
%       eval_id theta_1 theta_2 ... theta_<theta_len>
%   and this function returns it as a struct with fields eval_id, theta and
%   signature.
%
%   ok is false when the file is missing, unreadable, empty, unparseable or
%   carries the wrong number of values, so the caller can poll without a
%   try/catch. A short read that lands between the driver's write and its
%   flush therefore fails the length check and is retried on the next poll
%   rather than evaluated as a truncated theta.
%
%   eval_id is what pairs a result with the request that asked for it. Matching
%   on the theta values alone cannot do that: the acquisition maximiser can
%   propose the same point twice, and a float that has been through text is not
%   reliably equal to the one that was sent.
%
%   signature combines the file modification time with the content of the line.
%   It changes only when a new request arrives, which lets the caller skip a
%   file it has already served.

    req = struct("eval_id", NaN, "theta", [], "signature", "");
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
    if numel(vals) ~= theta_len + 1
        return
    end

    eval_id = vals(1);
    if ~isfinite(eval_id) || eval_id < 1 || eval_id ~= round(eval_id)
        return
    end

    req.eval_id = eval_id;
    req.theta = vals(2:end);

    d = dir(theta_txt);
    if isempty(d)
        req.signature = "missing|" + last_line;
    else
        req.signature = string(d.datenum) + "|" + last_line;
    end
    ok = true;
end
