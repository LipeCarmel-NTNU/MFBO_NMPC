function header = results_csv_header(theta_len)
%RESULTS_CSV_HEADER The one definition of the results CSV column order.
%
%   init_results_csv writes this header. check_results_header compares an
%   existing file against it. append_results_row writes its values in this
%   order. One definition is what makes the check work. A schema change that
%   misses one of the three sites fails at the next start. Without the check it
%   would show up later as a column that holds the wrong quantity.

    header = "eval_id,timestamp,phase,phi_vintage,z," + ...
        "SSE_measured,SSdU_measured,phi_SSE,phi_SSdU," + ...
        "SSE,SSdU,J,runtime_s,n_flag_not_one,phi_floored," + ...
        "wall_total_s,wall_cases_s,wall_phi_s,wall_build_s,wall_save_s";

    for k = 1:theta_len
        header = header + sprintf(",theta_%d", k);
    end
end
