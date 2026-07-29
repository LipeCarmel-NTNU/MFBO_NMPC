function append_results_row(results_csv, eval_id, ts, phase, out, theta)
%APPEND_RESULTS_ROW Append one evaluation to the results CSV.
%
%   The row is assembled from the out struct returned by simulate_nmpc, so the
%   measured costs, the divisors and the surrogate vintage travel together and
%   cannot fall out of step with one another.
%
%   Values are written with %.17g so that a double survives the round trip
%   through text. The line is built in memory and written by a single fprintf,
%   because the driver polls this file from another process: a short buffered
%   append lands whole, whereas a sequence of small writes can expose a partial
%   row to a reader that opens the file between them.

    line = sprintf("%d,%s,%s,%.17g,%.17g", ...
        eval_id, ts, phase, ...
        field_or_nan(out, "beta_vintage"), ...
        cfg_fidelity(out));

    line = line + sprintf(",%.17g,%.17g,%.17g,%.17g", ...
        field_or_nan(out, "SSE_measured"), field_or_nan(out, "SSdU_measured"), ...
        field_or_nan(out, "frac_SSE"),     field_or_nan(out, "frac_SSdU"));

    line = line + sprintf(",%.17g,%.17g,%.17g,%.17g,%.17g,%d", ...
        field_or_nan(out, "SSE"),  field_or_nan(out, "SSdU"), ...
        field_or_nan(out, "J"),    field_or_nan(out, "runtime_s"), ...
        field_or_nan(out, "n_flag_not_one"), ...
        logical_or_zero(out, "frac_floored"));

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
%FIELD_OR_NAN Read a scalar field, returning NaN when it is absent or empty.
% An absent field is written as NaN rather than skipped, so every row carries
% the same number of columns whatever the phase produced it.
    if isfield(s, name) && ~isempty(s.(name))
        v = double(s.(name));
    else
        v = NaN;
    end
end


function v = logical_or_zero(s, name)
%LOGICAL_OR_ZERO Read a flag field as 0 or 1, defaulting to 0.
    if isfield(s, name) && ~isempty(s.(name)) && s.(name)
        v = 1;
    else
        v = 0;
    end
end


function z = cfg_fidelity(out)
%CFG_FIDELITY The simulated fidelity, read back from the decoded configuration.
% Taking it from out.cfg rather than from theta(1) records the value the
% simulation actually ran at, after decode_theta clamped it into [0, 1].
    if isfield(out, "cfg") && isfield(out.cfg, "f")
        z = double(out.cfg.f);
    else
        z = NaN;
    end
end
