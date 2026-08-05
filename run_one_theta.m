function out = run_one_theta(theta_or_id, opts)
%RUN_ONE_THETA Rerun a single theta through the pipeline simulator and plot it.
%
%   out = run_one_theta(theta)        theta is a 1x12 vector, same layout as the
%                                     campaign: [f, theta_p, theta_m,
%                                     log10(Q), log10(Ru), log10(Rdu)].
%   out = run_one_theta(id)           id is a stored evaluation id
%                                     (out_<n1>_<n2>.mat -> n1n2); its theta is
%                                     loaded from the matching .mat.
%   out = run_one_theta()             uses benchmark_theta().
%
%   This calls the same simulate_nmpc the driver uses, so the result matches a
%   real evaluation. rng(1) is set first to reproduce the campaign noise draw.
%
%   Name-value options:
%     horizon       "full" (default, whole base.tf) or "fidelity" (tf = f*base.tf,
%                   what the BO loop actually evaluates)
%     extrapolate   scale partial costs by phi(z) (default false); needs a phi file
%     phi_coeffs    phi_coeffs.mat path (default results/surrogate/phi_coeffs.mat)
%     terminal_cost "lqr" (default), "zero" or "none"
%     x0            initial-condition cases, one row each (default: simulate_nmpc's)
%     sigma_y       measurement-noise std (default [0.001 0.1 0.1]; [0 0 0] = noiseless)
%     results_root  where to search for the id .mat (default "results")
%     plot          draw trajectories (default true)
%
%   Examples:
%     out = run_one_theta(20260802022215);              % rerun a stored eval, full horizon
%     out = run_one_theta(th, horizon="fidelity", extrapolate=true);  % reproduce a BO eval
%     out = run_one_theta([1 0 1 -1.6 0.65 -3 3 3 1 3 3 -3], sigma_y=[0 0 0]);

    arguments
        theta_or_id double = []
        opts.horizon (1,1) string {mustBeMember(opts.horizon, ["full" "fidelity"])} = "full"
        opts.extrapolate (1,1) logical = false
        opts.phi_coeffs (1,1) string = ""
        opts.terminal_cost (1,1) string {mustBeMember(opts.terminal_cost, ["lqr" "zero" "none"])} = "lqr"
        opts.x0 (:,:) double = []
        opts.sigma_y (1,:) double = [0.001 0.1 0.1]
        opts.results_root (1,1) string = "results"
        opts.plot (1,1) logical = true
    end

    %% Put the project on the path
    here = fileparts(mfilename("fullpath"));
    if ~isfolder(fullfile(here, "dependencies"))
        here = fileparts(here);   % allow the file to live one level down
    end
    addpath(genpath(here));

    %% Resolve theta (vector as-is, scalar id -> load, empty -> benchmark)
    if isempty(theta_or_id)
        theta = benchmark_theta();
        fprintf("Using benchmark_theta().\n");
    elseif isscalar(theta_or_id)
        [theta, src] = theta_from_id(theta_or_id, fullfile(here, opts.results_root));
        fprintf("Loaded theta of id %d from %s\n", theta_or_id, src);
    else
        theta = theta_or_id(:)';
    end
    if numel(theta) ~= 12
        error("run_one_theta:theta", "theta must have 12 elements, got %d.", numel(theta));
    end

    %% Build the shared base and run
    phi_path = opts.phi_coeffs;
    if opts.extrapolate && strlength(phi_path) == 0
        phi_path = fullfile(here, opts.results_root, "surrogate", "phi_coeffs.mat");
    end

    rng(1);
    base = nmpc_base(sigma_y = opts.sigma_y, phi_coeffs_path = phi_path);

    sim_args = {"horizon", opts.horizon, "extrapolate", opts.extrapolate, ...
                "terminal_cost", opts.terminal_cost};
    if ~isempty(opts.x0)
        sim_args = [sim_args, {"x0", opts.x0}];
    end
    out = simulate_nmpc(base, theta, sim_args{:});

    %% Report
    fprintf("\ntheta = [%s]\n", strtrim(sprintf("%.4g ", theta)));
    fprintf("horizon=%s  extrapolate=%d  terminal_cost=%s\n", ...
        opts.horizon, opts.extrapolate, opts.terminal_cost);
    fprintf("SSE=%.6g  SSdU=%.6g  J=%.6g\n", out.SSE, out.SSdU, out.J);
    fprintf("SSE_measured=%.6g  SSdU_measured=%.6g  n_flag_not_one=%d\n", ...
        out.SSE_measured, out.SSdU_measured, out.n_flag_not_one);

    if opts.plot
        plot_cases(out);
    end
end


function [theta, src] = theta_from_id(id, results_root)
%THETA_FROM_ID Find out_<n1>_<n2>.mat whose concatenated digits equal id.
    files = dir(fullfile(results_root, "**", "out_*.mat"));
    for k = 1:numel(files)
        tok = regexp(files(k).name, "^out_(\d+)_(\d+)\.mat$", "tokens", "once");
        if isempty(tok); continue; end
        if str2double([tok{1}, tok{2}]) == id
            src = fullfile(files(k).folder, files(k).name);
            S = load(src, "out");
            theta = S.out.theta(:)';
            return
        end
    end
    error("run_one_theta:id", "No out_*.mat with id %d under %s.", id, results_root);
end


function plot_cases(out)
%PLOT_CASES One figure per scenario: three states with setpoints, and inputs.
    labels = ["V", "X", "S"];
    for k = 1:numel(out.case)
        c = out.case(k);
        n = size(c.Y, 1);
        t = (0:n-1) * c.dt;
        figure("Name", sprintf("scenario %d", k));
        for s = 1:3
            subplot(4, 1, s);
            plot(t, c.Y(:, s), "b-", "LineWidth", 2, "DisplayName", "plant"); hold on;
            plot(t, c.Ysp(:, s), "r--", "LineWidth", 2, "DisplayName", "setpoint");
            ylabel(labels(s)); grid on; box on;
            if s == 1
                title(sprintf("scenario %d", k), "Interpreter", "none");
                legend("Location", "best");
            end
        end
        subplot(4, 1, 4);
        plot(t, c.U, "LineWidth", 2); grid on; box on;
        ylabel("u"); xlabel("time (h)");
    end
end
