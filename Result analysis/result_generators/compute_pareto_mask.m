function isPareto = compute_pareto_mask(J_track, J_TV)
%COMPUTE_PARETO_MASK Mark non-dominated tradeoff points (lower is better).
n = numel(J_track);
isPareto = true(n, 1);
for i = 1:n
    dominated = (J_track <= J_track(i)) & (J_TV <= J_TV(i)) & ...
                ((J_track < J_track(i)) | (J_TV < J_TV(i)));
    dominated(i) = false;
    if any(dominated)
        isPareto(i) = false;
    end
end
end
