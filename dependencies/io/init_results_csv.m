function init_results_csv(results_csv, theta_len)
%INIT_RESULTS_CSV Write the results CSV header if the file does not exist.
%
%   Columns, in order:
%     eval_id         monotonic evaluation index assigned by the driver
%     timestamp       run identifier, also the name of the out_<ts>.mat file
%     phase           DOE or OPT
%     beta_vintage    surrogate fit that scaled this row, NaN when unscaled
%     z               simulated fidelity, theta(1)
%     SSE_measured    tracking cost measured at z, before any scaling
%     SSdU_measured   control-variation cost measured at z, before any scaling
%     frac_SSE        divisor applied to SSE_measured
%     frac_SSdU       divisor applied to SSdU_measured
%     SSE, SSdU       full-horizon estimates passed to the optimiser
%     J               SSE + 1e4 * SSdU
%     runtime_s       wall-clock cost of the evaluation
%     n_flag_not_one  steps whose fmincon exit flag was not 1
%     frac_floored    1 when the 0.01 floor on either divisor was active
%     theta_1..theta_<theta_len>
%
%   The measured costs and the divisors sit next to the scaled ones so that any
%   row can be rescaled under a later surrogate fit without rerunning the
%   simulation. An existing file is left untouched so that a run can append to
%   it after an interruption; check_results_header verifies that the columns of
%   an existing file are the ones this version writes.

    if isfile(results_csv)
        return
    end

    ensure_parent_dir(results_csv);

    fid = fopen(results_csv, "w");
    if fid < 0
        error("Could not open results file for writing: %s", results_csv);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "%s\n", results_csv_header(theta_len));
end
