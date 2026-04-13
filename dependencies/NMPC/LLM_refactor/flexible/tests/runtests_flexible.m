function results = runtests_flexible()
    % Runner for the flexible-version test suite. Adds the flexible folder
    % and its subfolders (for the component classes) plus the tests folder
    % to the MATLAB path, runs all tests, then restores the path.

    here     = fileparts(mfilename('fullpath'));
    flex_dir = fileparts(here);

    added = { ...
        flex_dir, ...
        fullfile(flex_dir, 'model'), ...
        fullfile(flex_dir, 'integrator'), ...
        fullfile(flex_dir, 'constraints'), ...
        fullfile(flex_dir, 'cost'), ...
        fullfile(flex_dir, 'solver'), ...
        here};
    old = path;
    cleanup = onCleanup(@() path(old));
    for k = 1 : numel(added)
        addpath(added{k});
    end

    results = runtests(here);
    disp(results);
end
