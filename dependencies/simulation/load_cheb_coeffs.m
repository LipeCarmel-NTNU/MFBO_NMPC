function coeffs = load_cheb_coeffs(path)
%LOAD_CHEB_COEFFS Load the fidelity surrogate coefficients fitted on the
%initialisation runs.
%
%   Returns a struct with fields SSE and SSdU, each a column vector of
%   Chebyshev coefficients c_0..c_n for the cost fraction curve
%   frac(f) = sum_k c_k T_k(2*f - 1).
%
%   The file is written by
%   J surrogate/runtime_surrogate/fit_runtime_surrogate.py after the
%   initialisation phase. A missing file is an error: running BO against
%   coefficients from an earlier vintage would make the objective values
%   untraceable to the initialisation set that produced them.

    if ~isfile(path)
        error(['Fidelity surrogate coefficients not found: %s\n' ...
               'Run the initialisation phase (main_initialization.m), then fit the\n' ...
               'surrogate with J surrogate/runtime_surrogate/fit_runtime_surrogate.py\n' ...
               'before starting BO.'], path);
    end

    S = load(path);
    if ~isfield(S, "cheb_c_SSE") || ~isfield(S, "cheb_c_SSdU")
        error("%s must contain cheb_c_SSE and cheb_c_SSdU.", path);
    end

    coeffs = struct();
    coeffs.SSE = double(S.cheb_c_SSE(:));
    coeffs.SSdU = double(S.cheb_c_SSdU(:));

    if isempty(coeffs.SSE) || isempty(coeffs.SSdU)
        error("Coefficient vectors in %s are empty.", path);
    end
    if any(~isfinite(coeffs.SSE)) || any(~isfinite(coeffs.SSdU))
        error("Coefficient vectors in %s contain non-finite values.", path);
    end
end
