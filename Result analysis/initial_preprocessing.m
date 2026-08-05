
%% initial_preprocessing
% Build a flat, human-friendly view of the MFBO-NMPC results and save it to
% Result analysis/storage/preprocessed.mat.
%
% Terminology
%   case     : one full campaign (DOE + BO), e.g. case1
%   phase    : DOE or BO within a case
%   scenario : one setpoint scenario inside an evaluation (2 per evaluation)
%
% Output variables in preprocessed.mat
%   evals     : table, one row per BO evaluation (the object you filter and plot)
%   doe       : table, one row per DOE initialisation evaluation (same schema)
%   surrogate : table, one row per phi(z) vintage (coeffs, lambda, refit time)
%   checks    : table, per-BO-evaluation phi diagnostics for quick sanity checks
%   meta      : provenance of this preprocessing run
%
% Trajectories are not stored; retrieve them later with get_simulation(id).
%
% Run parse_registry.py first so that storage/surrogate_vintages.csv exists.
% The list of case subfolders lives in cases.txt (shared with the parser).

%% Paths and cases
here        = fileparts(mfilename('fullpath'));
repo_root   = fileparts(here);
results_root = fullfile(repo_root, 'results');
storage_dir = fullfile(here, 'storage');
if ~exist(storage_dir, 'dir'); mkdir(storage_dir); end

% Case subfolders to ingest, read from cases.txt (shared with parse_registry.py).
cases = read_cases(fullfile(here, 'cases.txt'));

%% Ingest evaluations
rows = table();
skipped   = strings(0, 1);
mismatches = strings(0, 1);

for ic = 1:numel(cases)
    case_name = cases{ic};
    bo_dir  = fullfile(results_root, case_name);
    doe_dir = fullfile(bo_dir, 'init');

    phase_dirs = {doe_dir, 'DOE'; bo_dir, 'BO'};
    for ip = 1:size(phase_dirs, 1)
        d     = phase_dirs{ip, 1};
        phase = phase_dirs{ip, 2};
        files = dir(fullfile(d, 'out_*.mat'));
        for k = 1:numel(files)
            fpath = fullfile(files(k).folder, files(k).name);
            try
                r = read_eval(fpath, case_name, phase);
            catch ME
                skipped(end+1, 1) = string(fpath) + " : " + string(ME.message); %#ok<AGROW>
                continue
            end
            rows = [rows; r]; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    error('initial_preprocessing:noData', 'No out_*.mat files found under %s', results_root);
end

%% Execution order (iter) per case and phase, then categoricals
rows.case  = categorical(rows.case);
rows.phase = categorical(rows.phase, {'DOE', 'BO'});
rows.iter  = zeros(height(rows), 1);
grp = findgroups(rows.case, rows.phase);   % iter counts within each case+phase
for g = 1:max(grp)
    idx = find(grp == g);
    [~, ord] = sort(rows.timestamp(idx));  % ascending execution order
    iter = zeros(numel(idx), 1);
    iter(ord) = (1:numel(idx))';
    rows.iter(idx) = iter;
end

%% Split into BO evaluations (evals) and DOE initialisation (doe), same schema
cols = { ...
    'id', 'case', 'iter', 'timestamp', ...
    'z', 'Np', 'Nc', 'Q1', 'Q2', 'Q3', ...
    'Ru1', 'Ru2', 'Ru3', 'Rdu1', 'Rdu2', 'Rdu3', ...
    'SSE', 'SSdU', 'SSE_measured', 'SSdU_measured', 'N', ...
    'n_flag0', 'n_flag1', 'n_flag2', 'n_flag_neg', ...
    't_total', 't_nmpc'};
evals = sortrows(rows(rows.phase == 'BO',  cols), {'case', 'iter'});
doe   = sortrows(rows(rows.phase == 'DOE', cols), {'case', 'iter'});

%% Per-evaluation phi diagnostics for the BO runs
% Expect any(phi_floored)==0 and any(~extrapolated)==0. DOE runs carry no phi.
checks = rows(rows.phase == 'BO', {'id', 'case', 'phi_vintage', 'phi_floored', 'extrapolated'});
checks = sortrows(checks, {'case', 'id'});

%% Surrogate vintages (produced by parse_registry.py)
csv_path = fullfile(storage_dir, 'surrogate_vintages.csv');
if exist(csv_path, 'file')
    surrogate = readtable(csv_path);
    if ismember('case', surrogate.Properties.VariableNames)
        surrogate.case = categorical(surrogate.case);
    end
else
    warning('initial_preprocessing:noSurrogate', ...
        'Missing %s. Run parse_registry.py first. surrogate left empty.', csv_path);
    surrogate = table();
end

%% Provenance
meta = struct();
meta.generated_at = datetime('now');
meta.repo_root    = repo_root;
meta.results_root = results_root;
meta.cases        = cases;
meta.n_evals      = height(evals);   % BO evaluations
meta.n_doe        = height(doe);     % DOE initialisation evaluations
meta.n_vintages   = height(surrogate);
meta.n_skipped    = numel(skipped);
meta.skipped      = skipped;
meta.mismatches   = mismatches;
meta.note         = ['Trajectories not stored: use get_simulation(id). ', ...
    'n_flag_neg counts EXITFLAG<0 summed over both scenarios; ', ...
    'if EXITFLAG is stored as an unsigned type this may read 0.'];

%% Save
out_path = fullfile(storage_dir, 'preprocessed.mat');
save(out_path, 'evals', 'doe', 'surrogate', 'checks', 'meta', '-v7.3');
fprintf('Saved %d BO evaluations, %d DOE, %d vintages -> %s\n', ...
    height(evals), height(doe), height(surrogate), out_path);
if ~isempty(skipped)
    fprintf('Skipped %d file(s); see meta.skipped\n', numel(skipped));
end


function cases = read_cases(path)
%READ_CASES Read case subfolder names from cases.txt (one per line, # comments).
if ~exist(path, 'file')
    error('read_cases:missing', 'Case list not found: %s', path);
end
lines = readlines(path);
lines = strip(lines);
keep  = lines ~= "" & ~startsWith(lines, "#");
cases = cellstr(lines(keep));
if isempty(cases)
    error('read_cases:empty', 'No cases listed in %s', path);
end
end


function r = read_eval(fpath, case_name, phase)
%READ_EVAL Extract one evaluation row from an out_*.mat file.
S   = load(fpath);
out = S.out;
cfg = out.cfg;

% id from filename out_<n1>_<n2>.mat -> concatenated integer
[~, name] = fileparts(fpath);
tok = regexp(name, '^out_(\d+)_(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('read_eval:badName', 'Unexpected filename: %s', name);
end
id        = int64(str2double([tok{1}, tok{2}]));
timestamp = datetime([tok{1}, tok{2}], 'InputFormat', 'yyyyMMddHHmmss');

% Decision variables: z from theta(1); horizons and weights from out.cfg
theta = out.theta(:)';
z  = theta(1);
Np = double(cfg.p);            % prediction horizon
Nc = double(cfg.m);            % control horizon
Qd  = diag(cfg.Q);
Rud = diag(cfg.Ru);
Rdud = diag(cfg.Rdu);

% Flag counts by EXITFLAG value, summed over scenarios
n0 = 0; n1 = 0; n2 = 0; nneg = 0;
for s = 1:numel(out.case)
    ef = double(out.case(s).EXITFLAG(:));
    n0   = n0   + sum(ef == 0);
    n1   = n1   + sum(ef == 1);
    n2   = n2   + sum(ef == 2);
    nneg = nneg + sum(ef < 0);
end

r = table( ...
    id, string(case_name), string(phase), timestamp, ...
    z, Np, Nc, Qd(1), Qd(2), Qd(3), ...
    Rud(1), Rud(2), Rud(3), Rdud(1), Rdud(2), Rdud(3), ...
    out.SSE, out.SSdU, out.SSE_measured, out.SSdU_measured, double(out.N), ...
    n0, n1, n2, nneg, ...
    out.wall_s.total, out.runtime_s, ...
    double(out.phi_vintage), logical(out.phi_floored), logical(out.extrapolated), ...
    'VariableNames', { ...
    'id', 'case', 'phase', 'timestamp', ...
    'z', 'Np', 'Nc', 'Q1', 'Q2', 'Q3', ...
    'Ru1', 'Ru2', 'Ru3', 'Rdu1', 'Rdu2', 'Rdu3', ...
    'SSE', 'SSdU', 'SSE_measured', 'SSdU_measured', 'N', ...
    'n_flag0', 'n_flag1', 'n_flag2', 'n_flag_neg', ...
    't_total', 't_nmpc', ...
    'phi_vintage', 'phi_floored', 'extrapolated'});
end
