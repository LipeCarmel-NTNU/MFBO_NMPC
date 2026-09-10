function gen_num_runtime_consistency(ctx)
%GEN_NUM_RUNTIME_CONSISTENCY Timing consistency figures as numbers.
%   Result: results/numerical results/runtime_consistency.txt
%   Per case: medians and sums of t_total and t_nmpc, the gap distribution,
%   Spearman correlations against failed solves, and a list of evaluations
%   whose timing is inconsistent (large positive gap, or negative from a
%   checkpoint resume).
T = ra_require(ctx, "timeline");
A = T.A; doeCount = T.doeCount; nCases = T.nCases; caseLabels = T.caseLabels;
numericalDir = ctx.numericalDir;

%% Numerical summary (console + txt)
txtPath = fullfile(numericalDir, 'runtime_consistency.txt');
fid = fopen(txtPath, 'w');
if fid == -1
    warning('runtime_consistency:txt', 'Could not open %s; console only.', txtPath);
    fid = 1;
end
emit(fid, 'Runtime consistency: t_total (wall_s.total) vs t_nmpc (runtime_s)');
emit(fid, 'gap = t_total - t_nmpc. Expected: small positive (minutes).');
emit(fid, 'gap >> 0 -> non-solver wall overhead; gap < 0 -> checkpoint resume artifact.');
emit(fid, '');

for k = 1:nCases
    T = A{k};
    emit(fid, '=== %s (%d evaluations: %d DOE + %d BO) ===', ...
        caseLabels(k), height(T), doeCount(k), height(T) - doeCount(k));
    emit(fid, '  t_total : median %.2f h | max %.2f h | sum %.1f h', ...
        median(T.totalH, 'omitnan'), max(T.totalH), sum(T.totalH, 'omitnan'));
    emit(fid, '  t_nmpc  : median %.2f h | max %.2f h | sum %.1f h', ...
        median(T.nmpcH, 'omitnan'), max(T.nmpcH), sum(T.nmpcH, 'omitnan'));
    emit(fid, '  gap     : median %.2f min | p95 %.2f min | max %.2f min | min %.2f min', ...
        median(T.gapMin, 'omitnan'), prctile_safe(T.gapMin, 95), ...
        max(T.gapMin), min(T.gapMin));
    emit(fid, '  gap/t_total: median %.2f %% | max %.2f %%', ...
        100 * median(T.gapMin * 60 ./ double(T.t_total), 'omitnan'), ...
        100 * max(T.gapMin * 60 ./ double(T.t_total)));

    % Correlations: is the wall time jump explained by solving or by the gap?
    rt = rank_corr(T.nBad, T.totalH);
    rn = rank_corr(T.nBad, T.nmpcH);
    rg = rank_corr(T.nBad, T.gapMin);
    emit(fid, '  Spearman rho vs failed solves: t_total %.3f | t_nmpc %.3f | gap %.3f', ...
        rt, rn, rg);

    % Offenders: gap above max(10 min, 10 %% of t_total), or negative gap.
    thr = max(10, 0.10 * double(T.t_total) / 60);
    bad = T.gapMin > thr | T.gapMin < -1;
    if any(bad)
        emit(fid, '  %d evaluation(s) with inconsistent timing:', nnz(bad));
        B = sortrows(T(bad, :), 'gapMin', 'descend');
        for i = 1:height(B)
            emit(fid, ['    id=%d | %s iter %3d | z=%.3f | t_total %6.2f h | ' ...
                't_nmpc %6.2f h | gap %8.2f min | failed solves %d%s'], ...
                B.id(i), string(B.phase_label(i)), B.iter(i), B.z(i), ...
                B.totalH(i), B.nmpcH(i), B.gapMin(i), B.nBad(i), ...
                resume_note(B.gapMin(i)));
        end
    else
        emit(fid, '  No inconsistent evaluations: t_total is a faithful runtime.');
    end
    emit(fid, '');
end

if fid ~= 1
    fclose(fid);
    fprintf('Wrote %s\n', txtPath);
end
end

function emit(fid, fmt, varargin)
%EMIT Print one line to the file and mirror it on the console.
line = sprintf(fmt, varargin{:});
fprintf(fid, '%s\n', line);
if fid ~= 1
    fprintf('%s\n', line);
end
end


function r = rank_corr(x, y)
%RANK_CORR Spearman rank correlation without toolbox dependencies.
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 3 || all(x == x(1)) || all(y == y(1))
    r = NaN;
    return
end
r = corrcoef(tied_rank(x), tied_rank(y));
r = r(1, 2);
end


function rk = tied_rank(v)
%TIED_RANK Average ranks with ties.
[~, ord] = sort(v);
rk = zeros(size(v));
rk(ord) = 1:numel(v);
u = unique(v);
for i = 1:numel(u)
    m = v == u(i);
    rk(m) = mean(rk(m));
end
end


function p = prctile_safe(v, q)
%PRCTILE_SAFE Percentile without the Statistics Toolbox.
v = sort(v(isfinite(v)));
if isempty(v)
    p = NaN;
    return
end
idx = max(1, min(numel(v), round(q / 100 * numel(v))));
p = v(idx);
end


function s = resume_note(gapMin)
%RESUME_NOTE Tag negative gaps as checkpoint-resume artifacts.
if gapMin < -1
    s = '  <- resumed run: t_total undercounts';
else
    s = '';
end
end
