function y = beta_eval(z, a, b)
%BETA_EVAL Fidelity cost fraction f(z) = I_z(a, b).
%
%   I_z(a, b) is the regularised incomplete beta function, increasing in z with
%   f(0) = 0 and f(1) = 1 for every positive a and b. The endpoints hold
%   exactly, so a full-fidelity evaluation is divided by exactly one and needs
%   no special case in the caller.
%
%   The shape parameters are fitted by
%   J surrogate/runtime_surrogate/fit_beta_surrogate.py. Note the argument
%   order: MATLAB's betainc takes (x, a, b) whereas scipy.special.betainc takes
%   (a, b, x). The two agree once the order is matched.
%
%   Arguments are clipped to the domain betainc accepts. The limits on (a, b)
%   are far wider than the fit bounds, so they act only on a corrupted
%   coefficient file, which load_beta_coeffs already rejects.

    ab_min = 1e-6;
    ab_max = 1e4;

    a = min(max(a, ab_min), ab_max);
    b = min(max(b, ab_min), ab_max);
    y = betainc(min(max(z, 0), 1), a, b);
end
