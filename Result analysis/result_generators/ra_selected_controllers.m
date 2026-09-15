function Sel = ra_selected_controllers(ctx, nMF, opts)
%RA_SELECTED_CONTROLLERS The controller set the selected_* figures draw.
%
%   Sel = ra_selected_controllers(ctx) returns a struct array, one element per
%   controller, ordered by INCREASING J_track:
%       .label     "BO_1", "MF_1", "RDU_50", ...
%       .arm       "BO", "MF" or "RDU"
%       .caseName  the campaign folder the run came from
%       .file      the .mat it was loaded from
%       .out       the loaded out struct
%       .z .SSE .SSdU .Np .Nc
%
%   Selection rule, data-driven so the set follows the campaigns instead of a
%   pinned list of timestamps:
%     - the nBO (default 1) lowest-J_track Pareto-optimal BO points of the
%       single-fidelity arm,
%     - the nMF (default 3) lowest-J_track BO points of the multi-fidelity arm,
%     - every damped-Rdu blend in results/rdu_damping/ (see ra_rdu_controllers).
%   DOE rows are excluded, exactly as they are from the frontier itself.
%
%   The single-fidelity Pareto set has two members; nBO = 1 keeps BO_1 and
%   drops BO_2, whose slot the Rdu blend now occupies. Pass nBO = Inf to get
%   the whole single-fidelity frontier back.
%
%   The multi-fidelity runs stopped at z < 1, so their out.T is SHORTER than
%   the full 10 h. Nothing here pads them; the figures draw each trace to its
%   own end and the truncation is meant to be visible.

arguments
    ctx struct
    nMF (1,1) double = 3
    opts.nBO (1,1) double = 1
    opts.rdu (1,1) logical = true
    opts.rduPattern (1,1) string = "rdu_*.mat"
end

F      = ra_require(ctx, "frontier");
names  = string(F.caseNames);
isBase = contains(lower(names), "baseline");
kBase  = find(isBase, 1);
kMF    = find(~isBase, 1);
if isempty(kBase) || isempty(kMF)
    error('ra_selected_controllers:arms', ...
        'Need one baseline arm and one multi-fidelity arm; cases.txt gave: %s', ...
        strjoin(names, ', '));
end

B = sortrows(F.Tp{kBase}, 'SSE');            % Tp is already the Pareto subset
B = B(1:min(opts.nBO, height(B)), :);
M = sortrows(F.E{kMF},    'SSE');
M = M(1:min(nMF, height(M)), :);

rows  = [B; M];
arm   = [repmat("BO", height(B), 1); repmat("MF", height(M), 1)];
cname = [repmat(names(kBase), height(B), 1); repmat(names(kMF), height(M), 1)];
label = arm + "_" + string([(1:height(B))'; (1:height(M))']);

Sel = repmat(ra_selected_template(), height(rows), 1);

for i = 1:height(rows)
    d = rows.timestamp(i);
    d.Format = 'yyyyMMdd_HHmmss';
    fpath = fullfile(ctx.resultsRoot, char(cname(i)), "out_" + string(d) + ".mat");
    if ~isfile(fpath)
        error('ra_selected_controllers:missing', ...
            '%s selected but its trends file is absent: %s', label(i), fpath);
    end
    L = load(fpath, 'out');

    Sel(i).label    = label(i);
    Sel(i).arm      = arm(i);
    Sel(i).caseName = cname(i);
    Sel(i).file     = string(fpath);
    Sel(i).out      = L.out;
    Sel(i).z        = double(rows.z(i));
    Sel(i).SSE      = double(rows.SSE(i));
    Sel(i).SSdU     = double(rows.SSdU(i));
    Sel(i).Np       = double(rows.Np(i));
    Sel(i).Nc       = double(rows.Nc(i));
end

if opts.rdu
    R = ra_rdu_controllers(ctx, opts.rduPattern);
    Sel = [Sel(:); reshape(R, [], 1)];
end

[~, ord] = sort([Sel.SSE]);                  % columns run by increasing J_track
Sel = Sel(ord);

for i = 1:numel(Sel)
    fprintf('%-7s %-14s z = %.4f  N = %4d  J_track = %10.2f  J_TV = %9.5f  Np = %2d  Nc = %d\n', ...
        Sel(i).label, Sel(i).caseName, Sel(i).z, double(Sel(i).out.N), ...
        Sel(i).SSE, Sel(i).SSdU, round(Sel(i).Np), round(Sel(i).Nc));
end
end
