function write_runtime_and_parameter_summary(E, Tp, caseLabels, ...
    runSummaryTables, combinedSummaryTable, outDir)
%WRITE_RUNTIME_AND_PARAMETER_SUMMARY Persist manuscript-facing runtime/tuning summaries.
nCases  = numel(E);
outPath = fullfile(outDir, 'runtime_and_params.txt');
fid = fopen(outPath, 'w');
if fid == -1
    warning('Unable to write numerical summary: %s', outPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for k = 1:nCases
    fprintf(fid, 'Runtime phase summary - %s:\n', char(caseLabels(k)));
    write_runtime_table(fid, runSummaryTables{k});
    fprintf(fid, '\n');
end
fprintf(fid, 'Runtime phase summary - All cases (combined):\n');
write_runtime_table(fid, combinedSummaryTable);
fprintf(fid, '\n');

% Mean fidelity z during the BO phase.
Eall  = vertcat(E{:});
meanZ = cellfun(@(T) mean(double(T.z), 'omitnan'), E);
fprintf('Mean fidelity z during BO evaluations:\n');
fprintf(fid, 'Mean fidelity z during BO evaluations:\n');
for k = 1:nCases
    fprintf('  %s: %.6g\n', caseLabels(k), meanZ(k));
    fprintf(fid, '  %s: %.6g\n', char(caseLabels(k)), meanZ(k));
end
meanZAll = mean(double(Eall.z), 'omitnan');
fprintf('  All cases (combined): %.6g\n', meanZAll);
fprintf(fid, '  All cases (combined): %.6g\n\n', meanZAll);

% Share of BO points with N_c = 1.
fprintf('Percentage of BO points with N_c = 1:\n');
fprintf(fid, 'Percentage of BO points with N_c = 1:\n');
nc1Lines = strings(nCases + 1, 1);
for k = 1:nCases
    n1 = nnz(double(E{k}.Nc) == 1);
    nT = height(E{k});
    nc1Lines(k) = sprintf('  %s: %d/%d (%.6g%%)', caseLabels(k), n1, nT, 100 * n1 / max(nT, 1));
    fprintf('%s\n', nc1Lines(k));
    fprintf(fid, '%s\n', nc1Lines(k));
end
n1 = nnz(double(Eall.Nc) == 1);
nT = height(Eall);
nc1Lines(end) = sprintf('  All cases (combined): %d/%d (%.6g%%)', n1, nT, 100 * n1 / max(nT, 1));
fprintf('%s\n', nc1Lines(end));
fprintf(fid, '%s\n\n', nc1Lines(end));

% Compact standalone txt for manuscript bookkeeping.
nc1Path = fullfile(outDir, 'optimization_nc1_share.txt');
fidNc1 = fopen(nc1Path, 'w');
if fidNc1 ~= -1
    cleanupNc1 = onCleanup(@() fclose(fidNc1)); %#ok<NASGU>
    fprintf(fidNc1, 'Percentage of BO points with N_c = 1:\n');
    fprintf(fidNc1, '%s\n', nc1Lines);
else
    warning('Unable to write N_c=1 share summary: %s', nc1Path);
end

% Per-case Pareto horizon composition.
for k = 1:nCases
    T = Tp{k};
    nRows   = height(T);
    countNc1 = nnz(double(T.Nc) == 1);
    countNp1 = nnz(double(T.Np) == 1);
    fprintf(fid, '%s Pareto counts (BO phase only, DOE excluded):\n', char(caseLabels(k)));
    fprintf(fid, '  N_c = 1: %d/%d\n', countNc1, nRows);
    fprintf(fid, '  N_p = 1: %d/%d\n', countNp1, nRows);

    ssdU    = double(T.SSdU);
    minSSdU = min(ssdU, [], 'omitnan');
    topMask = abs(ssdU - minSSdU) <= 10 * eps(max(1, abs(minSSdU)));
    fprintf(fid, '  Top (lowest SSdU) rows: %d\n', nnz(topMask));
    fprintf(fid, '  Top (lowest SSdU) with N_p = 1 and N_c = 1: %d\n\n', ...
        nnz(topMask & double(T.Np) == 1 & double(T.Nc) == 1));
end
fprintf('Saved: %s\n', outPath);
end
