function summaryTbl = display_runtime_phase_summary(Ddoe, Ebo, runName)
%DISPLAY_RUNTIME_PHASE_SUMMARY Runtime share spent in DOE vs BO (t_total).
doeMin   = sum(double(Ddoe.t_total), 'omitnan') / 60;
boMin    = sum(double(Ebo.t_total),  'omitnan') / 60;
totalMin = doeMin + boMin;

if totalMin > 0
    doePct = 100 * doeMin / totalMin;
    boPct  = 100 * boMin  / totalMin;
else
    doePct = NaN;
    boPct  = NaN;
end

phase             = {'DOE'; 'Optimisation'; 'Total'};
iterations        = [height(Ddoe); height(Ebo); height(Ddoe) + height(Ebo)];
runtime_min       = [doeMin; boMin; totalMin];
runtime_pct_total = [doePct; boPct; 100];

summaryTbl = table(phase, iterations, runtime_min, runtime_pct_total, ...
    'VariableNames', {'phase', 'iterations', 'runtime_min', 'runtime_pct_total'});

disp("Runtime phase summary - " + string(runName) + ":");
disp(summaryTbl);
end
