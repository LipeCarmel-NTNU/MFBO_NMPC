classdef test_smoke < matlab.unittest.TestCase
    % End-to-end smoke checks for the flexible composition-based NMPC.

    methods (TestClassSetup)
        function require_fmincon(tc)
            if isempty(which('fmincon'))
                tc.assumeFail('fmincon not available; skipping solve tests.');
            end
        end
    end

    methods (Test)
        %% Layout pack/unpack roundtrip
        function layout_roundtrip(tc)
            L = DecisionLayout();
            L.add('x', 5, 3);
            L.add('u', 2, 4);
            L.freeze();

            parts.x = rand(5, 3);
            parts.u = rand(2, 4);
            w = L.pack(parts, []);
            tc.verifyEqual(numel(w), 5*3 + 2*4);
            out = L.unpack(w, []);
            tc.verifyEqual(out.x, parts.x, 'AbsTol', 0);
            tc.verifyEqual(out.u, parts.u, 'AbsTol', 0);
        end

        function layout_with_scaling(tc)
            L = DecisionLayout();
            L.add('x', 3, 2);
            L.freeze();
            S = Scaling();
            S.set('x', [10 0.1]);

            parts.x = [1 2; 3 4; 5 6];
            w   = L.pack(parts, S);
            out = L.unpack(w, S);
            tc.verifyEqual(out.x, parts.x, 'AbsTol', 1e-12);
        end

        %% Build + single solve
        function single_solve_succeeds(tc)
            nmpc = build_basic_nmpc();
            [uk, x, u, info] = nmpc.solve([0 0], [0 0]);

            tc.verifyEqual(size(uk), [1 nmpc.nu]);
            tc.verifyEqual(size(x),  [nmpc.p + 1, nmpc.nx]);
            tc.verifyEqual(size(u),  [nmpc.m,     nmpc.nu]);
            tc.verifyGreaterThanOrEqual(info.flag, 0);
            tc.verifyTrue(isfinite(info.fval));
        end

        %% Closed-loop drives to setpoint
        function closed_loop_converges(tc)
            % Steady-state inputs that hold y_sp = [1 2] with f = diag(-1,-2)*x + I*u
            % are u_ss = [1 4]. Use them as u_sp so R doesn't bias.
            nmpc = build_basic_nmpc(u_sp=[1 4]);
            x = [0 0];
            u_prev = [0 0];
            for k = 1 : 80
                uk = nmpc.solve(x, u_prev);
                x  = nmpc.step(x, uk);
                u_prev = uk;
            end
            tc.verifyEqual(x, [1 2], 'AbsTol', 5e-2);
        end

        %% Input bounds respected
        function input_bounds_respected(tc)
            nmpc = build_basic_nmpc(umin=[-0.5 -0.5], umax=[0.5 0.5]);
            uk = nmpc.solve([0 0], [0 0]);
            tc.verifyLessThanOrEqual(uk, [0.5 0.5] + 1e-8);
            tc.verifyGreaterThanOrEqual(uk, [-0.5 -0.5] - 1e-8);
        end

        %% Scaling does not change continuity feasibility
        function scaling_preserves_continuity(tc)
            for sx = {[1 1], [10 0.1]}
                nmpc = build_basic_nmpc(x_scale=sx{1});
                [~, ~, ~, info] = nmpc.solve([0.5 -0.3], [0.1 0.2]);
                tc.verifyGreaterThanOrEqual(info.flag, 0);
            end
        end

        %% Terminal cost adds when provided
        function terminal_cost_opt_in(tc)
            nmpc = build_basic_nmpc(P=eye(2));
            [~, ~, ~, info] = nmpc.solve([0 0], [0 0]);
            tc.verifyGreaterThanOrEqual(info.flag, 0);
        end

        %% Soft state bounds: build and solve
        function soft_bounds_solve(tc)
            nmpc = build_basic_nmpc( ...
                soft_mask=[true false], rho_L1=10, ...
                Ymin=[-0.1 -5], Ymax=[0.1 5]);
            [~, x, ~, info] = nmpc.solve([0 0], [0 0]);
            tc.verifyGreaterThanOrEqual(info.flag, 0);
            % Hard state must stay within bounds.
            tc.verifyLessThanOrEqual(max(x(:,2)),  5 + 1e-8);
            tc.verifyGreaterThanOrEqual(min(x(:,2)), -5 - 1e-8);
        end

        %% Δu bound respected
        function dumax_respected(tc)
            nmpc = build_basic_nmpc(dumax=[0.05 0.05]);
            uk = nmpc.solve([0 0], [0 0]);
            % u_prev = 0, so |uk| ≤ dumax.
            tc.verifyLessThanOrEqual(abs(uk), [0.05 0.05] + 1e-8);
        end

        %% Warm start populates latest_wopt
        function warm_start_populated(tc)
            nmpc = build_basic_nmpc();
            tc.verifyTrue(isempty(nmpc.latest_wopt));
            nmpc.solve([0 0], [0 0]);
            tc.verifyFalse(isempty(nmpc.latest_wopt));
        end

        %% Logging
        function log_appends(tc)
            nmpc = build_basic_nmpc(log_enabled=true);
            nmpc.solve([0 0],   [0 0]);
            nmpc.solve([0.1 0.1], [0 0]);
            log = nmpc.get_log();
            tc.verifyEqual(numel(log), 2);
            tc.verifyTrue(isfield(log, 'w_opt'));
            tc.verifyTrue(isfield(log, 'flag'));
        end

        %% Builder validation
        function builder_requires_model(tc)
            tc.verifyError(@() NMPCBuilder().withHorizon(5, 2).build(), ...
                'NMPCBuilder:model');
        end

        function builder_requires_cost(tc)
            f = @(x, u) (-x(:) + u(:)).';
            b = NMPCBuilder() ...
                .withModel(ContinuousModel(f, 2, 2, RK4(0.1))) ...
                .withHorizon(5, 2);
            tc.verifyError(@() b.build(), 'NMPCBuilder:cost');
        end

        function builder_rejects_m_gt_p(tc)
            f = @(x, u) (-x(:) + u(:)).';
            b = NMPCBuilder() ...
                .withModel(ContinuousModel(f, 2, 2, RK4(0.1))) ...
                .withHorizon(3, 5) ...
                .addCost(QuadraticTracking(eye(2), eye(2), [0 0], [0 0]));
            tc.verifyError(@() b.build(), 'NMPCBuilder:horizon');
        end
    end
end
