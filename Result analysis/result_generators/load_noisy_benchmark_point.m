function [Jtrack, JTV, found] = load_noisy_benchmark_point(projectRoot)
%LOAD_NOISY_BENCHMARK_POINT Load noisy benchmark aggregate objectives (guarded).
benchPath = fullfile(projectRoot, 'results', 'benchmark_reference_controller', ...
    'benchmark_full_f1_same_noise_fix', 'out_benchmark.mat');
Jtrack = nan;
JTV = nan;
found = false;
if ~isfile(benchPath)
    warning('Noisy benchmark file not found; benchmark overlays skipped: %s', benchPath);
    return
end
S = load(benchPath, 'out');
if ~isfield(S, 'out') || ~isfield(S.out, 'SSE') || ~isfield(S.out, 'SSdU')
    warning('Noisy benchmark file missing required out.SSE/out.SSdU fields: %s', benchPath);
    return
end
Jtrack = double(S.out.SSE);
JTV = double(S.out.SSdU);
found = isfinite(Jtrack) && isfinite(JTV) && Jtrack > 0 && JTV > 0;
if ~found
    warning('Noisy benchmark values are invalid; benchmark overlays skipped: %s', benchPath);
end
end
