function write_runtime_table(fid, summaryTbl)
%WRITE_RUNTIME_TABLE Write one runtime summary table to disk.
for i = 1:height(summaryTbl)
    fprintf(fid, '  %s | iterations=%d | runtime_min=%.6g | runtime_pct_total=%.6g\n', ...
        char(string(summaryTbl.phase{i})), summaryTbl.iterations(i), ...
        summaryTbl.runtime_min(i), summaryTbl.runtime_pct_total(i));
end
end

%% ===================== PLOTTING HELPERS =====================
