function y = cheb_eval(x, c)
%CHEB_EVAL Evaluate a Chebyshev expansion sum_k c_k T_k(x).
%
%   c holds the coefficients c_0 to c_n, so the order follows from
%   numel(c) - 1. x is expected in [-1,1]; for the fidelity surrogate the
%   caller maps f in [0,1] with x = 2*f - 1.
%
%   The recurrence T_0 = 1, T_1 = x, T_{k+1} = 2*x*T_k - T_{k-1} matches
%   cheb_features in J surrogate/runtime_surrogate/fit_runtime_surrogate.py,
%   so coefficients fitted there can be evaluated here unchanged.

    c = c(:);
    if isempty(c)
        error("cheb_eval requires at least one coefficient.");
    end

    T_prev = ones(size(x));
    y = c(1) * T_prev;
    if numel(c) == 1
        return
    end

    T_curr = x;
    y = y + c(2) * T_curr;

    for k = 3:numel(c)
        T_next = 2 * x .* T_curr - T_prev;
        y = y + c(k) * T_next;
        T_prev = T_curr;
        T_curr = T_next;
    end
end
