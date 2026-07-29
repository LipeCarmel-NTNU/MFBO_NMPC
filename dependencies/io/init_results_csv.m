function init_results_csv(results_csv, theta_len)
%INIT_RESULTS_CSV Write the results CSV header if the file does not exist.
%
%   Columns, in order:
%     eval_id         evaluation index that the driver assigns
%     timestamp       run identifier, and the name of the out_<ts>.mat file
%     phase           DOE, OPT, FULL or BENCH
%     phi_vintage     surrogate fit that scaled this row, NaN when unscaled
%     z               simulated fidelity, theta(1)
%     SSE_measured    tracking cost measured at z, before any scaling
%     SSdU_measured   control-variation cost measured at z, before any scaling
%     phi_SSE         divisor applied to SSE_measured
%     phi_SSdU        divisor applied to SSdU_measured
%     SSE, SSdU       full-horizon estimates that the optimizer receives
%     J               SSE + 1e4 * SSdU
%     runtime_s       solver time summed over the control steps
%     n_flag_not_one  steps whose fmincon exit flag was not 1
%     phi_floored     1 when the 0.01 floor on either divisor acted
%     wall_total_s    wall time of the whole simulate_nmpc call
%     wall_cases_s    wall time of the case loop
%     wall_phi_s      wall time of the phi evaluation
%     wall_build_s    wall time of build_nmpc
%     wall_save_s     wall time of the .mat write
%     theta_1..theta_<theta_len>
%
%   The measured costs and the divisors sit next to the scaled costs. You can
%   therefore rescale any row under a later fit of phi without a new simulation.
%
%   The function leaves an existing file alone, so a run can append to it after
%   an interruption. check_results_header verifies that the columns of an
%   existing file are the columns that this version writes.

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
