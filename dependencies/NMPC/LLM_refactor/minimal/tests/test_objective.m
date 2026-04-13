classdef test_objective < matlab.unittest.TestCase
    % Verify each cost term contributes exactly the expected amount and that
    % optional terms are off by default.

    methods (Test)
        %% Tracking + input penalty only
        function tracking_and_input_only(tc)
            nmpc = build_basic_nmpc();
            % All states at setpoint, all inputs zero ⇒ cost = 0.
            x = repmat(nmpc.y_sp, nmpc.p + 1, 1);
            u = zeros(nmpc.m, nmpc.nu);
            ws = nmpc.pack_phys(x, u, []);
            J = nmpc.objfun(ws, zeros(1, nmpc.nu));
            tc.verifyEqual(J, 0, 'AbsTol', 1e-12);
        end

        function tracking_quadratic_in_error(tc)
            nmpc = build_basic_nmpc();
            x = zeros(nmpc.p + 1, nmpc.nx);   % off setpoint
            u = zeros(nmpc.m,     nmpc.nu);
            ws = nmpc.pack_phys(x, u, []);
            J = nmpc.objfun(ws, zeros(1, nmpc.nu));

            % Hand computation: p tracking errors of value -y_sp.
            e = -nmpc.y_sp.';
            J_expected = nmpc.p * (e.' * nmpc.Q * e);
            tc.verifyEqual(J, J_expected, 'AbsTol', 1e-10);
        end

        %% Terminal cost
        function terminal_cost_adds_when_P_set(tc)
            P = diag([3 7]);
            xspT = [4 -1];
            nmpc = build_basic_nmpc(P=P, x_sp_terminal=xspT);
            x = repmat(nmpc.y_sp, nmpc.p + 1, 1);   % zero tracking
            u = zeros(nmpc.m, nmpc.nu);
            ws = nmpc.pack_phys(x, u, []);
            J = nmpc.objfun(ws, zeros(1, nmpc.nu));

            e = (nmpc.y_sp - xspT).';
            tc.verifyEqual(J, e.' * P * e, 'AbsTol', 1e-10);
        end

        %% Δu penalty
        function du_cost_adds_when_S_set(tc)
            S = diag([2 4]);
            % Use R = 0 so we can isolate the Δu term.
            nmpc = build_basic_nmpc(S=S, R=zeros(2));
            x = repmat(nmpc.y_sp, nmpc.p + 1, 1);

            % Inputs constant and equal to u_prev ⇒ Δu = 0 throughout.
            u_prev = [1 2];
            u = repmat(u_prev, nmpc.m, 1);
            ws = nmpc.pack_phys(x, u, []);
            J = nmpc.objfun(ws, u_prev);
            tc.verifyEqual(J, 0, 'AbsTol', 1e-12);

            % Single step change at u(1), then constant: only Δu(1) ≠ 0.
            u_prev = [0 0];
            u = repmat([1 2], nmpc.m, 1);
            ws = nmpc.pack_phys(x, u, []);
            J = nmpc.objfun(ws, u_prev);
            duk = [1; 2];
            tc.verifyEqual(J, duk.' * S * duk, 'AbsTol', 1e-10);
        end

        %% Slack penalty (L1 + L2)
        function slack_penalty_L1(tc)
            rho = 7;
            nmpc = build_basic_nmpc(soft_mask=[true false], rho_L1=rho);
            x = repmat(nmpc.y_sp, nmpc.p + 1, 1);
            u = zeros(nmpc.m, nmpc.nu);
            s = ones(nmpc.p + 1, 1);
            ws = nmpc.pack_phys(x, u, s);
            J = nmpc.objfun(ws, zeros(1, nmpc.nu));
            tc.verifyEqual(J, rho * (nmpc.p + 1), 'AbsTol', 1e-10);
        end

        function slack_penalty_L2(tc)
            rho = 5;
            nmpc = build_basic_nmpc(soft_mask=[true false], rho_L2=rho);
            x = repmat(nmpc.y_sp, nmpc.p + 1, 1);
            u = zeros(nmpc.m, nmpc.nu);
            s = 2 * ones(nmpc.p + 1, 1);   % each slack squared = 4
            ws = nmpc.pack_phys(x, u, s);
            J = nmpc.objfun(ws, zeros(1, nmpc.nu));
            tc.verifyEqual(J, rho * 4 * (nmpc.p + 1), 'AbsTol', 1e-10);
        end
    end
end
