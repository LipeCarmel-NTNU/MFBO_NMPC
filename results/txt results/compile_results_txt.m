% COMPILE_RESULTS_TXT
% Combine all .txt files in this folder into:
%   "1 - Compiled results.txt"

clear; clc;

baseDir = fileparts(mfilename("fullpath"));
if strlength(string(baseDir)) == 0
    baseDir = pwd;
end

outName = "1 - Compiled results.txt";
outPath = fullfile(baseDir, outName);

txtFiles = dir(fullfile(baseDir, "*.txt"));
txtFiles = txtFiles(~strcmp({txtFiles.name}, outName));
[~, order] = sort(lower(string({txtFiles.name})));
txtFiles = txtFiles(order);

fidOut = fopen(outPath, "w");
if fidOut == -1
    error("Cannot open output file: %s", outPath);
end
cleanupObj = onCleanup(@() fclose(fidOut)); %#ok<NASGU>

for k = 1:numel(txtFiles)
    inPath = fullfile(txtFiles(k).folder, txtFiles(k).name);
    rawText = fileread(inPath);

    fprintf(fidOut, "===== %s =====\n", txtFiles(k).name);
    fprintf(fidOut, "%s", rawText);
    if ~endsWith(rawText, newline)
        fprintf(fidOut, "\n");
    end
    fprintf(fidOut, "\n");
end

fprintf("Compiled %d files into:\n%s\n", numel(txtFiles), outPath);
