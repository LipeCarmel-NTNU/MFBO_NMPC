function gen_tbl_pareto_candidates(ctx)
%GEN_TBL_PARETO_CANDIDATES The Pareto set of each case, printed as a table.
%   Result: console tables, one per case, sorted by descending J_track — the
%   candidate list the tuning-parameter table in the paper is drawn from.

    F = ra_require(ctx, "frontier");
    for k = 1:F.nCases
        fprintf('\n%s: %d BO evaluations, %d on the case frontier.\n', ...
            F.caseLabels(k), height(F.E{k}), height(F.Tp{k}));
        display_pareto_table(F.Tp{k}, F.caseLabels(k));
    end
end
