function y = phi_eval(z, a, b)
%PHI_EVAL Fidelity cost fraction phi(z) = I_z(a, b).
%
%   I_z(a, b) is the regularized incomplete beta function. phi increases with
%   z. phi(0) = 0 and phi(1) = 1 for every positive a and b. The end points hold
%   exactly, so a full-fidelity evaluation divides by exactly one. The caller
%   needs no special case at z = 1.
%
%   pipeline/phi_surrogate.py fits the shape parameters.
%
%   Note the argument order. MATLAB betainc takes (x, a, b).
%   scipy.special.betainc takes (a, b, x). The two agree after you match the
%   order.
%
%   This function clips its arguments to the domain that betainc accepts. The
%   limits on (a, b) are much wider than the fit bounds. They act only on a
%   corrupted coefficient file, which load_phi_coeffs already rejects.

    ab_min = 1e-6;
    ab_max = 1e4;

    a = min(max(a, ab_min), ab_max);
    b = min(max(b, ab_min), ab_max);
    y = betainc(min(max(z, 0), 1), a, b);
end
