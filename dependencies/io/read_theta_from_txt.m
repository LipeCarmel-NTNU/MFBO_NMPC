function [req, ok] = read_theta_from_txt(theta_txt, theta_len)
%READ_THETA_FROM_TXT Read one evaluation request from the inbox file.
%
%   The driver writes a single line
%       eval_id phase_code theta_1 theta_2 ... theta_<theta_len>
%   and this function returns it as a struct with the fields eval_id,
%   phase_code, theta and signature.
%
%   phase_code selects what the server does with the request:
%       0   design of experiments. Report the cost measured at the simulated
%           fidelity. Apply no surrogate.
%       1   optimization. Scale the measured cost by phi(z).
%   The driver sends the code with every request, so one server handles both
%   phases and no MATLAB script needs to know which phase a run is in.
%
%   ok is false when the file is missing, unreadable, empty, unparseable or
%   holds the wrong number of values. The caller can therefore poll without a
%   try/catch. A short read that lands between the write and the flush fails
%   the length check. The caller retries on the next poll instead of evaluating
%   a truncated theta.
%
%   eval_id pairs a result with the request that asked for it. Matching on the
%   theta values cannot do that. The acquisition maximizer can propose the same
%   point twice, and a float that has been through text is not reliably equal
%   to the one that was sent.
%
%   signature combines the file modification time with the content of the line.
%   It changes only when a new request arrives. The caller uses it to skip a
%   file that it has already served.

    req = struct("eval_id", NaN, "phase_code", NaN, "theta", [], "signature", "");
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
    if numel(vals) ~= theta_len + 2
        return
    end

    eval_id = vals(1);
    if ~isfinite(eval_id) || eval_id < 1 || eval_id ~= round(eval_id)
        return
    end

    phase_code = vals(2);
    if ~ismember(phase_code, [0 1])
        return
    end

    req.eval_id = eval_id;
    req.phase_code = phase_code;
    req.theta = vals(3:end);

    d = dir(theta_txt);
    if isempty(d)
        req.signature = "missing|" + last_line;
    else
        req.signature = string(d.datenum) + "|" + last_line;
    end
    ok = true;
end
