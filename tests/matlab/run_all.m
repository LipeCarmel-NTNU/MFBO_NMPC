function run_all()
%RUN_ALL Run every MATLAB test in this folder and report the result.
%
%   matlab -batch "run('tests/matlab/run_all.m')"
%
%   Each test_*.m is a function without arguments that prints one PASS line
%   per check and raises on the first failure. The tests run serially, without
%   a parallel pool, and take a few minutes. The Python tests are in
%   tests/python.

    here = fileparts(mfilename("fullpath"));
    addpath(here);
    files = dir(fullfile(here, "test_*.m"));
    failed = strings(0);
    for k = 1:numel(files)
        name = erase(files(k).name, ".m");
        fprintf("\n== %s\n", name);
        try
            feval(name);
        catch ME
            failed(end+1) = name; %#ok<AGROW>
            fprintf(2, "FAIL %s: %s\n", ME.identifier, ME.message);
        end
    end
    fprintf("\n%d test(s) failed", numel(failed));
    if isempty(failed)
        fprintf("\n");
    else
        fprintf(": %s\n", strjoin(failed, ", "));
        error("run_all:failed", "%d MATLAB test(s) failed.", numel(failed));
    end
end
