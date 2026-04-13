function results = runtests_minimal()
    % Run every unit test in this folder. Adds the parent (where NMPC.m lives)
    % and this folder to the path, then invokes runtests.
    %
    %   results = runtests_minimal();

    here   = fileparts(mfilename('fullpath'));
    parent = fileparts(here);

    addpath(parent, here);

    results = runtests(here, 'IncludeSubfolders', false);
    disp(results);
end
