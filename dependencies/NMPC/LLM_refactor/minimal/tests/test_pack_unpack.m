classdef test_pack_unpack < matlab.unittest.TestCase
    % Decision-vector layout, scaling roundtrip, bound construction.

    methods (Test)
        %% Roundtrip without slacks
        function roundtrip_no_slack(tc)
            nmpc = build_basic_nmpc();
            x = randn(nmpc.p + 1, nmpc.nx);
            u = randn(nmpc.m,     nmpc.nu);

            w = nmpc.pack_phys(x, u, []);
            [xb, ub, sb] = nmpc.unpack_phys(w);

            tc.verifyEqual(xb, x, 'AbsTol', 1e-12);
            tc.verifyEqual(ub, u, 'AbsTol', 1e-12);
            tc.verifyTrue(isempty(sb));
            tc.verifyEqual(numel(w), nmpc.len_x + nmpc.len_u);
        end

        %% Roundtrip with slacks
        function roundtrip_with_slack(tc)
            nmpc = build_basic_nmpc(soft_mask=[true false]);
            x = randn(nmpc.p + 1, nmpc.nx);
            u = randn(nmpc.m,     nmpc.nu);
            s = abs(randn(nmpc.p + 1, nmpc.n_soft));

            w = nmpc.pack_phys(x, u, s);
            [xb, ub, sb] = nmpc.unpack_phys(w);

            tc.verifyEqual(xb, x, 'AbsTol', 1e-12);
            tc.verifyEqual(ub, u, 'AbsTol', 1e-12);
            tc.verifyEqual(sb, s, 'AbsTol', 1e-12);
            tc.verifyEqual(numel(w), nmpc.len_x + nmpc.len_u + nmpc.len_s);
        end

        %% Scaling actually scales
        function scaling_changes_internal_w(tc)
            nmpc = build_basic_nmpc(x_scale=[2 5], u_scale=[3 7]);
            x = repmat([2 5], nmpc.p + 1, 1);     % equals x_scale
            u = repmat([3 7], nmpc.m,     1);     % equals u_scale
            w = nmpc.pack_phys(x, u, []);

            % Internally everything should be 1.
            tc.verifyEqual(w, ones(numel(w), 1), 'AbsTol', 1e-12);

            % Unpack returns physical values.
            [xb, ub] = nmpc.unpack_phys(w);
            tc.verifyEqual(xb, x, 'AbsTol', 1e-12);
            tc.verifyEqual(ub, u, 'AbsTol', 1e-12);
        end

        %% Bounds
        function bounds_have_right_dimensions(tc)
            nmpc = build_basic_nmpc(soft_mask=[true false]);
            [wL, wU] = nmpc.bounds_scaled();
            n = nmpc.len_x + nmpc.len_u + nmpc.len_s;
            tc.verifyEqual(numel(wL), n);
            tc.verifyEqual(numel(wU), n);
            % Slack lower bound is zero, upper is +Inf.
            tail = numel(wL) - nmpc.len_s + 1 : numel(wL);
            tc.verifyEqual(wL(tail), zeros(nmpc.len_s, 1));
            tc.verifyTrue(all(isinf(wU(tail))));
        end

        function soft_states_have_relaxed_box_bounds(tc)
            nmpc = build_basic_nmpc(soft_mask=[true false]);
            [wL, wU] = nmpc.bounds_scaled();
            xL = reshape(wL(1 : nmpc.len_x), [], nmpc.nx);
            xU = reshape(wU(1 : nmpc.len_x), [], nmpc.nx);

            % Soft column should be ±Inf; hard column should match Ymin/Ymax.
            tc.verifyTrue(all(isinf(xL(:, 1)) & xL(:, 1) < 0));
            tc.verifyTrue(all(isinf(xU(:, 1)) & xU(:, 1) > 0));
            tc.verifyEqual(xL(:, 2), repmat(nmpc.Ymin(2), nmpc.p + 1, 1));
            tc.verifyEqual(xU(:, 2), repmat(nmpc.Ymax(2), nmpc.p + 1, 1));
        end
    end
end
