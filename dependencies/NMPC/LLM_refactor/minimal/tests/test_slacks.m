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

        %% Tight Ymax forces non-zero slacks (otherwise infeasible)
        function slacks_activate_when_hard_bounds_infeasible(tc)
            % Setpoint y_sp = [1 2]; pin Ymax(1) BELOW it so reaching the
            % setpoint requires violating the hard upper bound.
            nmpc = build_basic_nmpc( ...
                Ymax      = [0.2  10], ...
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
                'Soft bound should be violated (slack > 0) once the closed loop drives x past Ymax.');
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
    end
end
