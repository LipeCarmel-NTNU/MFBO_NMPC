function gen_tbl_runtime_summary(ctx)
%GEN_TBL_RUNTIME_SUMMARY Phase runtime and frontier parameters per case.
%   Result: results/numerical results/ runtime and parameter summary, plus the
%   per-case and combined DOE-vs-BO phase tables on the console.

    F = ra_require(ctx, "frontier");
    P = ra_require(ctx, "preprocessed");

    tables = cell(F.nCases, 1);
    for k = 1:F.nCases
        tables{k} = display_runtime_phase_summary(F.D{k}, F.E{k}, F.caseLabels(k));
    end
    combined = display_runtime_phase_summary(P.doe, P.evals, "All cases (combined)");

    write_runtime_and_parameter_summary(F.E, F.Tp, F.caseLabels, ...
        tables, combined, ctx.numericalDir);
end
