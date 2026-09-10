function display_pareto_table(Tp, caseLabel)
%DISPLAY_PARETO_TABLE Print BO-phase Pareto controller settings for one case.
tuningCols = ["id", "iter", "SSE", "SSdU", "z", "Np", "Nc", "t_total", ...
    "Q1", "Q2", "Q3", "Ru1", "Ru2", "Ru3", "Rdu1", "Rdu2", "Rdu3"];
tuningCols = tuningCols(ismember(tuningCols, string(Tp.Properties.VariableNames)));
disp(caseLabel + " Pareto frontier with tuning weights (BO phase only, DOE excluded):");
disp(Tp(:, tuningCols));
end
