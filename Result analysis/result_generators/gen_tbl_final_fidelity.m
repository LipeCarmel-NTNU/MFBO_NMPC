function gen_tbl_final_fidelity(ctx)
%GEN_TBL_FINAL_FIDELITY Metrics for the frontier controllers re-run at z = 1.
%   Result: tables under results/numerical results/, settling time and IAE per
%   state for the Pareto candidates evaluated at full fidelity.
%   Guarded the same way as gen_fig_refined_frontier: needs the z=1
%   re-evaluation outputs on disk.

    F = ra_require(ctx, "frontier");
    F = ra_select_cases(F, "mf");   % nothing to refine in a run already at z = 1
    report_final_frontier_f1_metrics(F.E, F.caseNames, ctx.repo_root, ctx.numericalDir);
end
