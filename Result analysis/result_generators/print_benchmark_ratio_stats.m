function print_benchmark_ratio_stats(label, Jtrack, JTV, superiorMask, benchJtrack, benchJTV)
%PRINT_BENCHMARK_RATIO_STATS Print benchmark-relative ratios for superior controllers.
if nargin < 4 || isempty(superiorMask) || ~any(superiorMask)
    fprintf('    %s ratios: n/a (no strictly superior controllers)\n', label);
    return
end
JtrackSup = double(Jtrack(superiorMask));
JTVSup = double(JTV(superiorMask));
trackPct = 100 * mean(JtrackSup / benchJtrack, 'omitnan');
tvPct = 100 * mean(JTVSup / benchJTV, 'omitnan');
trackTimes = mean(benchJtrack ./ JtrackSup, 'omitnan');
tvTimes = mean(benchJTV ./ JTVSup, 'omitnan');
fprintf('    %s mean (J_better/J_bench): J_track=%.1f%%, J_TV=%.1f%%\n', label, trackPct, tvPct);
fprintf('    %s mean benchmark higher factor: J_track=%.1fx, J_TV=%.1fx\n', label, trackTimes, tvTimes);
end
