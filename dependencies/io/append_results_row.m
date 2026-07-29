function append_results_row(results_csv, eval_id, ts, phase, out, theta, wall_save_s)
%APPEND_RESULTS_ROW Append one evaluation to the results CSV.
%
%   The function builds the row from the out struct that simulate_nmpc returns.
%   The measured costs, the divisors and the surrogate vintage therefore travel
%   together and cannot fall out of step.
%
%   wall_save_s is the wall time of the .mat write. Pass it when the caller
%   measured it. The function writes NaN when the caller omits it.
%
%   The function writes each value with %.17g, so a double survives the round
%   trip through text. It builds the line in memory and writes it with one
%   fprintf. The driver polls this file from another process. One short
%   buffered append lands whole. A sequence of small writes can show a partial
%   row to a reader that opens the file between them.

    if nargin < 7
        wall_save_s = NaN;
    end

    line = sprintf("%d,%s,%s,%.17g,%.17g", ...
        eval_id, ts, phase, ...
        field_or_nan(out, "phi_vintage"), ...
        cfg_fidelity(out));

    line = line + sprintf(",%.17g,%.17g,%.17g,%.17g", ...
        field_or_nan(out, "SSE_measured"), field_or_nan(out, "SSdU_measured"), ...
        field_or_nan(out, "phi_SSE"),      field_or_nan(out, "phi_SSdU"));

    line = line + sprintf(",%.17g,%.17g,%.17g,%.17g,%.17g,%d", ...
        field_or_nan(out, "SSE"),  field_or_nan(out, "SSdU"), ...
        field_or_nan(out, "J"),    field_or_nan(out, "runtime_s"), ...
        field_or_nan(out, "n_flag_not_one"), ...
        logical_or_zero(out, "phi_floored"));

    line = line + sprintf(",%.17g,%.17g,%.17g,%.17g,%.17g", ...
        wall_or_nan(out, "total"), wall_or_nan(out, "cases"), ...
        wall_or_nan(out, "phi"),   wall_or_nan(out, "build_nmpc"), ...
        wall_save_s);

    theta = theta(:).';
    for k = 1:numel(theta)
        line = line + sprintf(",%.17g", theta(k));
    end

    fid = fopen(results_csv, "a");
    if fid < 0
        error("Could not open results file for appending: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "%s\n", line);
end


function v = field_or_nan(s, name)
%FIELD_OR_NAN Read a scalar field. Return NaN when it is absent or empty.
% The function writes NaN and does not skip the column. Every row therefore
% holds the same number of columns, whatever phase produced it.
    if isfield(s, name) && ~isempty(s.(name))
        v = double(s.(name));
    else
        v = NaN;
    end
end


function v = wall_or_nan(s, name)
%WALL_OR_NAN Read one stage time from out.wall_s. Return NaN when it is absent.
    if isfield(s, "wall_s") && isstruct(s.wall_s) && isfield(s.wall_s, name)
        v = double(s.wall_s.(name));
    else
        v = NaN;
    end
end


function v = logical_or_zero(s, name)
%LOGICAL_OR_ZERO Read a flag field as 0 or 1. Default to 0.
    if isfield(s, name) && ~isempty(s.(name)) && s.(name)
        v = 1;
    else
        v = 0;
    end
end


function z = cfg_fidelity(out)
%CFG_FIDELITY The simulated fidelity, read from the decoded configuration.
% The function takes the value from out.cfg and not from theta(1). It therefore
% records the fidelity that the simulation ran at, after decode_theta clamped it
% into [0, 1].
    if isfield(out, "cfg") && isfield(out.cfg, "f")
        z = double(out.cfg.f);
    else
        z = NaN;
    end
end
