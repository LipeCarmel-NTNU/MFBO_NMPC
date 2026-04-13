classdef test_solve < matlab.unittest.TestCase
    % End-to-end checks. Requires Optimization Toolbox (fmincon).

    methods (TestClassSetup)
        function require_fmincon(tc)
            if isempty(which('fmincon'))
                tc.assumeFail('fmincon not available; skipping solve tests.');
            end
        end
    end

    methods (Test)
        %% Single solve produces a feasible step
        function single_solve_succeeds(tc)
            nmpc = build_basic_nmpc();
            [uk, x, u, info] = nmpc.solve([0 0], [0 0]);

            tc.verifyEqual(size(uk), [1 nmpc.nu]);
            tc.verifyEqual(size(x),  [nmpc.p + 1, nmpc.nx]);
            tc.verifyEqual(size(u),  [nmpc.m,     nmpc.nu]);
            tc.verifyGreaterThanOrEqual(info.flag, 0);
            tc.verifyTrue(isfinite(info.fval));
        end

        %% Closed-loop drives state to setpoint
        function closed_loop_converges_to_setpoint(tc)
            % Steady-state inputs that hold y_sp = [1 2] are u_ss = [1 4]
            % (since A = diag(-1,-2), B = I). Use them as u_sp so the R
            % penalty doesn't bias the steady state.
            nmpc = build_basic_nmpc(u_sp=[1 4]);
            x = [0 0];
            u_prev = [0 0];
            n_steps = 80;
            for k = 1 : n_steps
                uk = nmpc.solve(x, u_prev);
                x  = nmpc.step(x, uk);
                u_prev = uk;
            end
            tc.verifyEqual(x, nmpc.y_sp, 'AbsTol', 5e-2, ...
                'Closed loop did not converge close to the setpoint.');
        end

        %% Bounds are honoured
        function input_bounds_respected(tc)
            % Tight umax should clip the optimal input.
            nmpc = build_basic_nmpc(umax=[0.5 0.5], umin=[-0.5 -0.5]);
            uk = nmpc.solve([0 0], [0 0]);
            tc.verifyLessThanOrEqual(uk, nmpc.umax + 1e-8);
            tc.verifyGreaterThanOrEqual(uk, nmpc.umin - 1e-8);
        end

        %% Warm start: latest_wopt populated after a successful solve
        function warm_start_state_is_populated(tc)
            nmpc = build_basic_nmpc();
            tc.verifyTrue(isempty(nmpc.latest_wopt));
            nmpc.solve([0 0], [0 0]);
            tc.verifyFalse(isempty(nmpc.latest_wopt));
            tc.verifyGreaterThanOrEqual(nmpc.latest_flag, 0);
        end

        %% Logging
        function log_appends_when_enabled(tc)
            nmpc = build_basic_nmpc(log_enabled=true);
            nmpc.solve([0 0], [0 0]);
            nmpc.solve([0.1 0.1], [0 0]);
            tc.verifyEqual(numel(nmpc.log), 2);
            tc.verifyTrue(isfield(nmpc.log, 'w_opt'));
            tc.verifyTrue(isfield(nmpc.log, 'flag'));
        end
    end
end
