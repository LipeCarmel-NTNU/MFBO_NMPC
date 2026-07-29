function header = results_csv_header(theta_len)
%RESULTS_CSV_HEADER The single definition of the results CSV column order.
%
%   init_results_csv writes it, check_results_header compares an existing file
%   against it, and append_results_row writes values in this order. Keeping one
%   definition is what makes the check meaningful: a schema change that misses
%   one of the three sites is caught on the next start rather than discovered
%   later in a column that silently holds the wrong quantity.

    header = "eval_id,timestamp,phase,beta_vintage,z," + ...
        "SSE_measured,SSdU_measured,frac_SSE,frac_SSdU," + ...
        "SSE,SSdU,J,runtime_s,n_flag_not_one,frac_floored";

    for k = 1:theta_len
        header = header + sprintf(",theta_%d", k);
    end
end
