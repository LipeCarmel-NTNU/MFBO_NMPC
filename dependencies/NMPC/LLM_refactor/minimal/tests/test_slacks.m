classdef test_slacks < matlab.unittest.TestCase
    % Soft state bounds via linear inequalities, slack penalty, and L1
    % exact-penalty recovery.

    methods (TestClassSetup)
        function require_fmincon(tc)
            if isempty(which('fmincon'))
                tc.assumeFail('fmincon not available; skipping slack tests.');
            end
        end
    end

    methods (Test)
        %% When the hard problem is feasible, slacks should be ~0 (L1 exact)
        function L1_exact_penalty_keeps_slacks_zero_when_feasible(tc)
            nmpc = build_basic_nmpc( ...
                soft_mask = [true false], ...
                rho_L1    = 1e4);
            nmpc.solve([0 0], [0 0]);

            [~, ~, s] = nmpc.unpack_phys(nmpc.latest_wopt);
            tc.verifyEqual(norm(s, Inf), 0, 'AbsTol', 1e-6, ...
                'Slacks should be zero when the hard problem is feasible.');
        end

        %% Tight Xmax forces non-zero slacks (otherwise infeasible)
        function slacks_activate_when_hard_bounds_infeasible(tc)
            % Setpoint x_sp = [1 2]; pin Xmax(1) BELOW it so reaching the
            % setpoint requires violating the hard upper bound.
            nmpc = build_basic_nmpc( ...
                Xmax      = [0.2  10], ...
                soft_mask = [true false], ...
                rho_L1    = 1);

            x = [0 0];
            u_prev = [0 0];
            for k = 1 : 30
                uk = nmpc.solve(x, u_prev);
                x  = nmpc.step(x, uk);
                u_prev = uk;
            end

            [~, ~, s] = nmpc.unpack_phys(nmpc.latest_wopt);
            tc.verifyGreaterThan(max(s(:)), 1e-3, ...
                'Soft bound should be violated (slack > 0) once the closed loop drives x past Xmax.');
        end

        %% Warm scalar solve: measured infeasibility is carried by slack
        function warm_start_infeasible_feedback_state_sets_slack_exactly(tc)
            xdot = @(x, u) -x(1) + u(1);
            nmpc = NMPC( ...
                xdot      = xdot, ...
                nx        = 1, ...
                nu        = 1, ...
                Ts        = 0.1, ...
                p         = 5, ...
                m         = 2, ...
                x_sp      = 0, ...
                u_sp      = 0, ...
                Q         = 1, ...
                R         = 0.1, ...
                Xmin      = -1, ...
                Xmax      = 1, ...
                umin      = -2, ...
                umax      = 2, ...
                soft_mask = true, ...
                rho_L1    = 1e3);

            nmpc.solve(0, 0);          % feasible point, stored as warm start
            [~, x, ~, info] = nmpc.solve(2, 0);
            [~, ~, s] = nmpc.unpack_phys(nmpc.latest_wopt);

            tc.verifyGreaterThanOrEqual(info.flag, 0);
            tc.verifyEqual(x(1), 2, 'AbsTol', 1e-8);
            tc.verifyEqual(s(1), 1, 'AbsTol', 1e-8, ...
                'Initial slack should equal x_init - Xmax exactly.');
        end

        %% Linear soft-bound rows have the right shape
        function soft_rows_dimensions(tc)
            nmpc = build_basic_nmpc(soft_mask=[true true]);
            [A, b] = nmpc.linear_ineq([0 0]);
            % 2 (upper+lower) × (p+1) × n_soft rows.
            expected_rows = 2 * (nmpc.p + 1) * nmpc.n_soft;
            tc.verifyEqual(size(A, 1), expected_rows);
            tc.verifyEqual(size(A, 2), nmpc.len_x + nmpc.len_u + nmpc.len_s);
            tc.verifyEqual(numel(b),   expected_rows);
            tc.verifyEqual(size(b), [expected_rows 1]);
        end

        function no_soft_no_linear_rows(tc)
            nmpc = build_basic_nmpc();        % no soft mask, no dumax
            [A, b] = nmpc.linear_ineq([0 0]);
            tc.verifyEqual(size(A, 1), 0);
            tc.verifyEqual(numel(b),   0);
        end

        %% dumax adds linear rows
        function dumax_rows_added(tc)
            nmpc = build_basic_nmpc(dumax=[0.1 0.1]);
            [A, b] = nmpc.linear_ineq([0 0]);
            % 2 directions × m × nu rows.
            expected = 2 * nmpc.m * nmpc.nu;
            tc.verifyEqual(size(A, 1), expected);
            tc.verifyEqual(numel(b),   expected);
        end

        function infinite_dumax_adds_no_rows(tc)
            nmpc = build_basic_nmpc(dumax=[Inf Inf]);
            [A, b] = nmpc.linear_ineq([0 0]);
            tc.verifyEqual(size(A, 1), 0);
            tc.verifyEqual(numel(b), 0);
        end

        function mixed_finite_and_infinite_dumax_only_adds_finite_rows(tc)
            nmpc = build_basic_nmpc(dumax=[0.1 Inf]);
            [A, b] = nmpc.linear_ineq([0 0]);
            expected = 2 * nmpc.m;
            tc.verifyEqual(size(A, 1), expected);
            tc.verifyEqual(numel(b), expected);
            tc.verifyTrue(all(isfinite(b)));
        end
    end
end
