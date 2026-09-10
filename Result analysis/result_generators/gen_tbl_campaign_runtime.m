function tbl = gen_tbl_campaign_runtime(ctx)
%GEN_TBL_CAMPAIGN_RUNTIME Cost of each campaign, as a table for the paper.
%
%   Result: results/numerical results/campaign_runtime.txt and .tex
%
%   One row per case: solver time t_nmpc summed over the DOE phase and over
%   the BO phase, the campaign total, and the mean cost of a BO evaluation.
%   t_nmpc is the comparable quantity across campaigns: it excludes checkpoint
%   I/O and, unlike t_total, it is accumulated across a resumed evaluation.
%   n is carried in the table because campaigns stopped on a wall-clock budget
%   have fewer BO evaluations than ones stopped on an iteration count, and the
%   per-evaluation column is the only fair comparison between them.

    F = ra_require(ctx, "frontier");

    caseLabel = strings(F.nCases, 1);
    nDoe = zeros(F.nCases, 1);  hDoe = zeros(F.nCases, 1);
    nBo  = zeros(F.nCases, 1);  hBo  = zeros(F.nCases, 1);
    for k = 1:F.nCases
        caseLabel(k) = F.caseLabels(k);
        nDoe(k) = height(F.D{k});
        nBo(k)  = height(F.E{k});
        hDoe(k) = sum(double(F.D{k}.t_nmpc), 'omitnan') / 3600;
        hBo(k)  = sum(double(F.E{k}.t_nmpc), 'omitnan') / 3600;
    end
    hTotal   = hDoe + hBo;
    hPerEval = hBo ./ max(nBo, 1);

    tbl = table(caseLabel, nDoe, hDoe, nBo, hBo, hTotal, hPerEval, ...
        'VariableNames', {'case', 'n_doe', 'doe_h', 'n_bo', 'bo_h', 'total_h', 'bo_h_per_eval'});
    disp(tbl);

    lines = strings(0, 1);
    lines(end+1) = "Campaign cost in NMPC solver time (t_nmpc = out.runtime_s).";
    lines(end+1) = "DOE and BO summed separately; per-evaluation column is the fair";
    lines(end+1) = "comparison when the BO budgets differ.";
    lines(end+1) = "";
    lines(end+1) = sprintf('%-28s %7s %9s %7s %9s %9s %11s', ...
        'campaign', 'n DOE', 'DOE (h)', 'n BO', 'BO (h)', 'total (h)', 'h/BO eval');
    for k = 1:F.nCases
        lines(end+1) = sprintf('%-28s %7d %9.2f %7d %9.2f %9.2f %11.3f', ...
            caseLabel(k), nDoe(k), hDoe(k), nBo(k), hBo(k), hTotal(k), hPerEval(k));
    end

    txtPath = fullfile(ctx.numericalDir, 'campaign_runtime.txt');
    write_lines(txtPath, lines);

    % LaTeX body only: no caption or label, so the paper keeps control of both.
    tex = strings(0, 1);
    tex(end+1) = '\begin{tabular}{lrrrrrr}';
    tex(end+1) = '\hline';
    tex(end+1) = 'Campaign & $n_{\mathrm{DOE}}$ & DOE (h) & $n_{\mathrm{BO}}$ & BO (h) & Total (h) & h/eval \\';
    tex(end+1) = '\hline';
    for k = 1:F.nCases
        tex(end+1) = sprintf('%s & %d & %.2f & %d & %.2f & %.2f & %.3f \\\\', ...
            latex_escape(caseLabel(k)), nDoe(k), hDoe(k), nBo(k), hBo(k), hTotal(k), hPerEval(k));
    end
    tex(end+1) = '\hline';
    tex(end+1) = '\end{tabular}';

    texPath = fullfile(ctx.numericalDir, 'campaign_runtime.tex');
    write_lines(texPath, tex);
end

function write_lines(path, lines)
%WRITE_LINES One string per line, with the format string kept free of content.
    fid = fopen(path, 'w');
    if fid == -1
        warning('gen_tbl_campaign_runtime:write', 'Could not open %s.', path);
        return
    end
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines(i));
    end
    fclose(fid);
    fprintf('Wrote %s\n', path);
end

function s = latex_escape(s)
%LATEX_ESCAPE Underscores in a case label would otherwise start a subscript.
    s = replace(string(s), "_", "\_");
end
