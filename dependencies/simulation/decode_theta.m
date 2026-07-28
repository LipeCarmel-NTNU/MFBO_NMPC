function cfg = decode_theta(theta, nx, nu)
%DECODE_THETA Decode a theta vector into an NMPC configuration.
%
%   theta = [f, theta_p, theta_m, q_exp(1:nx), ru_exp(1:nu), rdu_exp(1:nu)]
%   where q_exp, ru_exp and rdu_exp are log10 of the diagonal entries.
%
%   Horizons are encoded as
%       m = theta_m + 1
%       p = theta_p + m
%
%   f is the fidelity, clamped into [0,1] as a guard against invalid values
%   arriving from the external optimiser.

    theta = theta(:).';
    expected_len = 1 + 2 + nx + nu + nu;
    if numel(theta) ~= expected_len
        error("theta must have length %d.", expected_len)
    end

    k = 1;
    cfg = struct();

    cfg.f = theta(k); k = k + 1;

    theta_p = theta(k); k = k + 1;
    theta_m = theta(k); k = k + 1;

    cfg.m = theta_m + 1;
    cfg.p = theta_p + cfg.m;

    q_exp   = theta(k:k+nx-1); k = k + nx;
    ru_exp  = theta(k:k+nu-1); k = k + nu;
    rdu_exp = theta(k:k+nu-1);

    cfg.Q   = diag(10.^q_exp);
    cfg.Ru  = diag(10.^ru_exp);
    cfg.Rdu = diag(10.^rdu_exp);

    cfg.f = max(0, min(1, cfg.f));
end
