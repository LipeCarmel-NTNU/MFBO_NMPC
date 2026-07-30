function coeffs = load_phi_coeffs(path, opts)
%LOAD_PHI_COEFFS Load the fitted fidelity surrogate phi(z) = I_z(a, b).
%
%   The function returns a struct with the fields SSE and SSdU. Each field holds
%   the shape parameters a and b and the L2 strength lambda that produced them.
%   The struct also holds the vintage number and the file modification time.
%
%   Every results row keeps the vintage. You can therefore trace a stored
%   objective to the exact fit that scaled it.
%
%   pipeline/phi_surrogate.py writes the file. It
%   renames a temporary file over this path. A reader therefore sees the
%   previous vintage or the new one. The rename can still collide with a read on
%   Windows, so this function retries for a short time.
%
%   Name-value options:
%     required   raise an error when the file is absent (default true)
%     retries    number of attempts before the function gives up (default 5)
%     retry_s    pause between attempts, in seconds (default 0.2)

    arguments
        path (1,1) string
        opts.required (1,1) logical = true
        opts.retries (1,1) double {mustBePositive} = 5
        opts.retry_s (1,1) double {mustBeNonnegative} = 0.2
    end

    coeffs = [];

    if ~isfile(path)
        if opts.required
            error(['Fidelity surrogate coefficients not found: %s\n' ...
                   'Run the initialization phase (main_initialization.m). The\n' ...
                   'driver then fits phi and writes this file before the BO phase.'], path);
        end
        return
    end

    last_err = [];
    for attempt = 1:opts.retries
        try
            S = load(path);
            coeffs = validate_payload(S, path);
            return
        catch ME
            last_err = ME;
            pause(opts.retry_s);
        end
    end

    rethrow(last_err);
end


function coeffs = validate_payload(S, path)
%VALIDATE_PAYLOAD Check the loaded struct and build the coefficient record.
% A malformed file raises an error. The function applies no fallback. A run that
% continued against a partial or stale fit would change the meaning of every
% objective that follows, and no record would show it.

    targets = ["SSE", "SSdU"];
    coeffs = struct();

    for k = 1:numel(targets)
        t = targets(k);
        fa = "a_" + t;
        fb = "b_" + t;
        if ~isfield(S, fa) || ~isfield(S, fb)
            error("%s must contain %s and %s.", path, fa, fb);
        end

        a = double(S.(fa));
        b = double(S.(fb));
        if ~isscalar(a) || ~isscalar(b)
            error("%s: %s and %s must be scalars.", path, fa, fb);
        end
        if ~isfinite(a) || ~isfinite(b) || a <= 0 || b <= 0
            error("%s: %s = %g and %s = %g must be finite and positive.", ...
                path, fa, a, fb, b);
        end

        entry = struct("a", a, "b", b, "lambda", NaN);
        fl = "lambda_" + t;
        if isfield(S, fl)
            entry.lambda = double(S.(fl));
        end
        coeffs.(t) = entry;
    end

    if ~isfield(S, "vintage")
        error("%s must contain vintage.", path);
    end
    coeffs.vintage = double(S.vintage);
    if ~isscalar(coeffs.vintage) || ~isfinite(coeffs.vintage)
        error("%s: vintage must be a finite scalar.", path);
    end

    if isfield(S, "created_at")
        coeffs.created_at = string(S.created_at);
    else
        coeffs.created_at = "";
    end

    d = dir(path);
    if isempty(d)
        coeffs.file_datenum = NaN;
    else
        coeffs.file_datenum = d.datenum;
    end
    coeffs.path = string(path);
end
