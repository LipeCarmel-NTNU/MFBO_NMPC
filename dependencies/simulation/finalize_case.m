function case_out = finalize_case(base, x0, case_id, T, noise, Y, Y_meas, Ysp, U, RUNTIME, EXITFLAG, i_last)
%FINALIZE_CASE Assemble one case output and its cost sums.
%
%   Costs are accumulated over the first i_last samples, so a partially
%   completed case can be summarised at a checkpoint. State 1 carries a weight
%   of 10 in the tracking error, matching E = E.*[10 1 1] in the original
%   scripts.
%
%   partial_SSE and partial_SSdU are the per-step contributions. Their
%   cumulative sums are the cost trajectories that
%   J surrogate/runtime_surrogate/fit_runtime_surrogate.py fits the fidelity
%   surrogate to, so they are always stored, at every fidelity.

    Y_eval = Y(1:i_last, :);
    Ysp_eval = Ysp(1:i_last, :);
    U_eval = U(1:i_last, :);

    E = Y_eval - Ysp_eval;
    E = E .* [10 1 1];
    partial_SSE = sum(E.^2, 2);

    dU = diff(U_eval, 1, 1);
    partial_SSdU = sum(dU.^2, 2);

    case_out = struct();
    case_out.case_id = case_id;
    case_out.x0 = x0;

    case_out.Y = Y;
    case_out.Y_meas = Y_meas;
    case_out.Ysp = Ysp;
    case_out.U = U;

    case_out.noise = noise;
    case_out.dt = base.dt;
    case_out.tf = T(min(i_last, numel(T)));
    case_out.i_last = i_last;

    % Cost at the simulated fidelity. A caller that extrapolates to the full
    % horizon overwrites SSE and SSdU and keeps these as the measured values.
    case_out.SSE = sum(partial_SSE);
    case_out.SSdU = sum(partial_SSdU);
    case_out.partial_SSE = partial_SSE;
    case_out.partial_SSdU = partial_SSdU;

    case_out.RUNTIME = RUNTIME;
    case_out.EXITFLAG = EXITFLAG;
    case_out.runtime_s = sum(RUNTIME(1:i_last));

    % Steps whose fmincon exit flag was not 1. Non-empty means the optimiser
    % stopped for another reason and the corresponding actions deserve a look.
    flags = EXITFLAG(1:i_last);
    case_out.n_flag_not_one = sum(~(flags == 1));
    case_out.flag_steps = find(~(flags == 1));
end
