function coeffs = load_beta_coeffs(path, opts)
%LOAD_BETA_COEFFS Load the fitted fidelity surrogate f(z) = I_z(a, b).
%
%   Returns a struct with fields SSE and SSdU, each holding the shape
%   parameters a and b and the L2 strength lambda that produced them, together
%   with the vintage number and the file modification time. The vintage is
%   copied into every results row, so a stored objective can always be traced
%   to the exact fit that scaled it.
%
%   The file is written by
%   J surrogate/runtime_surrogate/fit_beta_surrogate.py, which publishes it by
%   renaming a temporary file over this path. A reader therefore sees either the
%   previous vintage or the new one. The rename can still collide with a
%   concurrent read on Windows, so a short bounded retry covers the transient.
%
%   Name-value options:
%     required   error when the file is missing (default true)
%     retries    attempts before giving up (default 5)
%     retry_s    pause between attempts in seconds (default 0.2)

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
                   'Run the initialisation phase (main_initialization.m), then fit\n' ...
                   'the surrogate with\n' ...
                   '  python "J surrogate/runtime_surrogate/fit_beta_surrogate.py" results/init\n' ...
                   'before starting BO.'], path);
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
%VALIDATE_PAYLOAD Check the loaded struct and assemble the coefficient record.
% A malformed file is an error rather than a fallback: continuing against a
% partially written or stale fit would silently change the meaning of every
% objective produced afterwards.

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
