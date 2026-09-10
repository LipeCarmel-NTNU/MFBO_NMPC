function main_results(varargin)
%MAIN_RESULTS Produce every result in the manuscript, one generator at a time.
%
%   main_results                       runs the whole chain
%   main_results('list')               prints the steps and stops
%   main_results('data')               only the preprocessing steps
%   main_results('tables')             only the numeric results
%   main_results('figures')            only the figures
%   main_results('runtime_cumulative') one step, prefix optional
%
%   Every step lives in result_generators/ and has one objective. A step is
%   either a partial preprocessing stage, which writes what it derived into
%   Result analysis/storage/ for the later steps to read, or a result, which
%   writes a figure into results/graphical_results/ or a table into
%   results/numerical results/. Nothing is recomputed twice: the per-case
%   split and the Pareto masks are derived once by gen_data_frontier, and the
%   execution timeline once by gen_data_timeline.
%
%   Prerequisite, outside MATLAB, whenever a campaign is added or extended:
%       python "Result analysis/parse_registry.py"
%   which refreshes cases.txt and storage/surrogate_vintages.csv. The case
%   list decides what appears in every figure, so a new campaign is a data
%   change and not a code change.
%
%   A step that raises is reported and the run continues, so one missing
%   input cannot cost you the rest of the results. The guarded steps
%   (refined frontier, final fidelity) skip themselves with a warning when
%   the z=1 re-evaluation outputs are not on disk.

    ctx   = ra_context();
    steps = local_steps();

    if nargin > 0 && any(strcmpi(varargin{1}, {'list', '-list', '--list'}))
        print_steps(steps);
        return
    end

    sel = select_steps(steps, varargin);
    fprintf('main_results: %d step(s)\n', numel(sel));

    failed = strings(0, 1);
    for i = sel
        fprintf('\n---- %-26s %s\n', steps(i).name, steps(i).about);
        t0 = tic;
        try
            steps(i).fn(ctx);
            fprintf('     done in %.1f s\n', toc(t0));
        catch ME
            failed(end+1, 1) = string(steps(i).name); %#ok<AGROW>
            fprintf(2, '     FAILED after %.1f s: %s\n', toc(t0), ME.message);
            fprintf(2, '     %s\n', getReport(ME, 'basic', 'hyperlinks', 'off'));
        end
    end

    fprintf('\nFigures  -> %s\n', ctx.graphicsDir);
    fprintf('Tables   -> %s\n', ctx.numericalDir);
    fprintf('Storage  -> %s\n', ctx.storageDir);
    if isempty(failed)
        fprintf('All %d step(s) completed.\n', numel(sel));
    else
        fprintf(2, '%d step(s) failed: %s\n', numel(failed), strjoin(failed, ', '));
    end
end


%% ===================== the step list =====================

function steps = local_steps()
%LOCAL_STEPS The one place that says what a full results run consists of.
    rows = {
    % name                      group      function                            objective
      'preprocess'              'data'     @run_preprocessing                  'ingest every out_*.mat -> storage/preprocessed.mat'
      'frontier'                'data'     @gen_data_frontier                  'per-case split and Pareto masks -> storage/frontier.mat'
      'timeline'                'data'     @gen_data_timeline                  'DOE-then-BO execution timeline -> storage/timeline.mat'
      'campaign_runtime'        'tables'   @gen_tbl_campaign_runtime           'cost of each campaign in t_nmpc, for the paper'
      'runtime_summary'         'tables'   @gen_tbl_runtime_summary            'phase runtimes and frontier parameters'
      'pareto_candidates'       'tables'   @gen_tbl_pareto_candidates          'the Pareto set of each case'
      'runtime_consistency_num' 'tables'   @gen_num_runtime_consistency        't_total vs t_nmpc, per case, as numbers'
      'final_fidelity'          'tables'   @gen_tbl_final_fidelity             'settling time and IAE at z = 1 (guarded)'
      'objective_space_z'       'figures'  @gen_fig_objective_space_z          'objective space per case, shaded by z'
      'runtime_iteration'       'figures'  @gen_fig_runtime_iteration          'per-evaluation runtime and z over iterations'
      'pareto_combined'         'figures'  @gen_fig_pareto_combined            'pooled samples and per-case frontiers'
      'runtime_cumulative'      'figures'  @gen_fig_runtime_cumulative         'cumulative runtime per case'
      'best_so_far'             'figures'  @gen_fig_best_so_far                'running best of each objective over the BO iterations'
      'refined_frontier'        'figures'  @gen_fig_refined_frontier           'low-fidelity frontier vs its z=1 refinement (guarded)'
      'runtime_consistency'     'figures'  @gen_fig_runtime_consistency        'solve time vs wall time, gap, failed solves'
      'surrogate_coefficients'  'figures'  @gen_fig_surrogate_coefficients     'fitted a, b, lambda per vintage'
      'surrogate_fit_quality'   'figures'  @gen_fig_surrogate_fit_quality      'fit loss and refit cost per vintage'
      'surrogate_phi_curves'    'figures'  @gen_fig_surrogate_phi_curves       'phi(z) per published vintage'
    };
    steps = struct('name', rows(:, 1), 'group', rows(:, 2), ...
                   'fn', rows(:, 3), 'about', rows(:, 4));
end


%% ===================== selection and reporting =====================

function idx = select_steps(steps, args)
%SELECT_STEPS Resolve the arguments to step indices, in list order.
    if isempty(args)
        idx = 1:numel(steps);
        return
    end
    names  = string({steps.name});
    groups = string({steps.group});
    keep   = false(1, numel(steps));
    for a = 1:numel(args)
        want = string(args{a});
        want = erase(want, ["gen_fig_", "gen_tbl_", "gen_num_", "gen_data_"]);
        if any(groups == want)
            keep = keep | (groups == want);
        elseif any(names == want)
            keep = keep | (names == want);
        else
            error('main_results:unknownStep', ...
                'No step or group named "%s". Run main_results(''list'').', want);
        end
    end
    idx = find(keep);
end

function print_steps(steps)
%PRINT_STEPS What a full run does, in order.
    fprintf('\nmain_results steps (run one with main_results(''<name>''))\n\n');
    lastGroup = '';
    for i = 1:numel(steps)
        if ~strcmp(steps(i).group, lastGroup)
            fprintf('  [%s]\n', steps(i).group);
            lastGroup = steps(i).group;
        end
        fprintf('    %-26s %s\n', steps(i).name, steps(i).about);
    end
    fprintf('\n');
end


%% ===================== the one step that is a script =====================

function run_preprocessing(ctx)
%RUN_PREPROCESSING Rebuild storage/preprocessed.mat from the results tree.
% initial_preprocessing is a script, so it runs here in a local workspace of
% its own rather than in the caller's. It reads cases.txt, so refresh that
% with parse_registry.py first if a campaign was added.
    casesFile = fullfile(ctx.here, 'cases.txt');
    if ~isfile(casesFile)
        error('main_results:noCases', ...
            'Missing %s. Run: python "%s"', casesFile, ...
            fullfile(ctx.here, 'parse_registry.py'));
    end
    initial_preprocessing;
end
