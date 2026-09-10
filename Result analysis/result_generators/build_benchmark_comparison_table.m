function Tcmp = build_benchmark_comparison_table(Jtrack, JTV, benchJtrack, benchJTV)
%BUILD_BENCHMARK_COMPARISON_TABLE Per-controller benchmark-relative metrics.
Jtrack = double(Jtrack(:));
JTV = double(JTV(:));
n = numel(Jtrack);

ratioTrackPct = nan(n, 1);
ratioTVPct = nan(n, 1);
timesTrack = nan(n, 1);
timesTV = nan(n, 1);
isSuperior = false(n, 1);

valid = isfinite(Jtrack) & isfinite(JTV) & (Jtrack > 0) & (JTV > 0);
ratioTrackPct(valid) = 100 * (Jtrack(valid) / benchJtrack);
ratioTVPct(valid) = 100 * (JTV(valid) / benchJTV);
timesTrack(valid) = benchJtrack ./ Jtrack(valid);
timesTV(valid) = benchJTV ./ JTV(valid);
isSuperior(valid) = (Jtrack(valid) < benchJtrack) & (JTV(valid) < benchJTV);

Tcmp = table(ratioTrackPct, ratioTVPct, timesTrack, timesTV, isSuperior, ...
    'VariableNames', {'J_track_ratio_pct', 'J_TV_ratio_pct', ...
    'J_track_bench_higher_x', 'J_TV_bench_higher_x', 'strictly_better'});
Tcmp = Tcmp(Tcmp.strictly_better, :);
Tcmp = Tcmp(:, {'J_track_ratio_pct', 'J_TV_ratio_pct', 'J_track_bench_higher_x', 'J_TV_bench_higher_x'});
if ~isempty(Tcmp)
    Tcmp = sortrows(Tcmp, 'J_TV_ratio_pct', 'ascend');
end
end

%% ===================== GUARDED: FINAL f=1 METRICS =====================
