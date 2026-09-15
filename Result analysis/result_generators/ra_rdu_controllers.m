function R = ra_rdu_controllers(ctx, pattern)
%RA_RDU_CONTROLLERS The damped-Rdu blends written by main_rdu_damping.
%
%   R = ra_rdu_controllers(ctx) scans results/rdu_damping/ for files matching
%   "rdu_*.mat" and returns one struct per BLEND run, in the field layout
%   ra_selected_controllers uses, so the two sets concatenate.
%
%   Each file holds `synth` (the synthesis result) and `runs` (the simulated
%   thetas). Only the BLEND run is taken; MF_1 and BO_1 are re-simulated there
%   only when SIM_MF / SIM_BO are switched on and they already have campaign
%   trends files of their own.
%
%   The label carries the loop-gain sacrifice: synth.EPS_X = 0.50 gives
%   "RDU_50". Several eps files therefore appear as several columns, ordered
%   with everything else by increasing J_track.
%
%   This reads results only. It writes nothing and moves nothing, so it is safe
%   to call while main_rdu_damping is still running another eps.

arguments
    ctx struct
    pattern (1,1) string = "rdu_*.mat"
end

R = ra_selected_template();
R(:) = [];                   % keep the field layout, drop the placeholder entry

rduDir = fullfile(ctx.resultsRoot, 'rdu_damping');
if ~isfolder(rduDir)
    return
end

files = dir(fullfile(rduDir, pattern));
files = files(~[files.isdir]);
if isempty(files)
    return
end

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    L = load(fpath, 'synth', 'runs');
    if ~isfield(L, 'runs') || isempty(L.runs)
        warning('ra_rdu_controllers:noRuns', 'No runs in %s; skipped.', fpath);
        continue
    end

    labels = string({L.runs.label});
    k = find(labels == "BLEND", 1);
    if isempty(k)
        warning('ra_rdu_controllers:noBlend', 'No BLEND run in %s; skipped.', fpath);
        continue
    end

    o   = L.runs(k).out;
    cfg = o.cfg;

    e = ra_selected_template();
    e.label    = sprintf('RDU_%d', round(100 * double(L.synth.EPS_X)));
    e.arm      = "RDU";
    e.caseName = "rdu_damping";
    e.file     = string(fpath);
    e.out      = o;
    e.z        = double(cfg.f);
    e.SSE      = double(o.SSE);
    e.SSdU     = double(o.SSdU);
    e.Np       = double(cfg.p);
    e.Nc       = double(cfg.m);

    R(end + 1) = e; %#ok<AGROW>
end

R = R(:);
end
