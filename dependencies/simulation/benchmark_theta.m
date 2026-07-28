function theta = benchmark_theta()
%BENCHMARK_THETA Theta of the benchmark reference controller.
%
%   Benchmark definition: m = 6, p = 61, Q = diag([10 1 1]), Ru = 0,
%   Rdu = diag([10 10 10]). The terminal cost P = 0 is not part of theta and
%   is set by the caller.
%
%   Layout: [f, theta_p, theta_m, log10(Q_diag), log10(Ru_diag), log10(Rdu_diag)]

    f = 1;
    m = 6;
    p = 61;
    theta_m = m - 1;
    theta_p = p - m;

    q_diag = [10 1 1];
    ru_exp = -1000 * ones(1, 3);    % 10^-1000 underflows to 0 in double precision
    rdu_diag = [10 10 10];

    theta = [f, theta_p, theta_m, log10(q_diag), ru_exp, log10(rdu_diag)];
end
