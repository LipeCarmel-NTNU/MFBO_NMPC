function lim = padded_log_limits(vals, padFactor)
%PADDED_LOG_LIMITS Data-driven limits for a positive log axis.
if nargin < 2 || ~isfinite(padFactor) || padFactor <= 1
    padFactor = 1.25;
end
vals = vals(isfinite(vals) & (vals > 0));
if isempty(vals)
    lim = [1e-2, 1e2];
    return
end
lim = [min(vals) / padFactor, max(vals) * padFactor];
end
