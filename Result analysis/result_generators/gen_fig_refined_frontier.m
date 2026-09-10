function gen_fig_refined_frontier(ctx)
%GEN_FIG_REFINED_FRONTIER Low-fidelity frontier against its z=1 refinement.
%   Result: results/graphical_results/refined_frontier_change.png/.pdf
%   Paper:  fig:refined
%
%   Guarded: skips with a warning unless the re-evaluation CSVs exist at
%   results/final_fidelity_same_noise/<case>_full_f1_same_noise/results_full.csv
%
%   The panel-b comparison has no meaning for a campaign whose evaluations
%   already ran at z = 1: there is nothing to refine.

    F = ra_require(ctx, "frontier");
    F = ra_select_cases(F, "mf");   % nothing to refine in a run already at z = 1
    run_refined_frontier_change(F.E, F.caseNames, F.caseLabels, ctx.repo_root, ...
        ctx.graphicsDir, ctx.plotColors, ctx.caseMarkers, ctx.accentColor, ctx.fontSize);
end
