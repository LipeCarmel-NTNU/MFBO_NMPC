function [nSuperior, superiorMask] = count_benchmark_strict_superior(Jtrack, JTV, benchJtrack, benchJTV)
%COUNT_BENCHMARK_STRICT_SUPERIOR Strictly better than benchmark in both objectives.
superiorMask = isfinite(Jtrack) & isfinite(JTV) & ...
    (Jtrack < benchJtrack) & (JTV < benchJTV);
nSuperior = nnz(superiorMask);
end
