function S = ra_require(ctx, name)
%RA_REQUIRE Load one artifact from Result analysis/storage, or say what to run.
%
%   S = ra_require(ctx, 'preprocessed') returns the evals/doe/surrogate tables.
%   S = ra_require(ctx, 'frontier')     returns the per-case split and masks.
%   S = ra_require(ctx, 'timeline')     returns the per-case execution timeline.
%
%   The generators read these instead of recomputing them, so a partial
%   preprocessing step runs once per campaign and not once per figure.

    path = fullfile(ctx.storageDir, name + ".mat");
    if ~isfile(path)
        switch name
            case "preprocessed", how = "parse_registry.py, then initial_preprocessing";
            case "frontier",     how = "gen_data_frontier(ctx)";
            case "timeline",     how = "gen_data_timeline(ctx)";
            otherwise,           how = "the generator that writes it";
        end
        error('ra_require:missing', 'Missing %s. Run %s first.', path, how);
    end
    S = load(path);
end
