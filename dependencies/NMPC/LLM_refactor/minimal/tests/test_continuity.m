classdef test_continuity < matlab.unittest.TestCase
    % confun should report ~zero continuity residual when w is built by
    % rolling out the model from x_init, and nonzero otherwise.

    methods (Test)
        function rolled_out_trajectory_satisfies_continuity(tc)
            nmpc = build_basic_nmpc();
            x_init = [0.5 -0.3];
            u_init = [0.1 0.2];

            % guess_from_initial returns a SCALED w (pack_phys handles the
            % conversion), which is exactly what confun expects.
            w0 = nmpc.guess_from_initial(x_init, u_init);
            [c, ceq] = nmpc.confun(w0, x_init);

            tc.verifyTrue(isempty(c));
            tc.verifyEqual(norm(ceq, Inf), 0, 'AbsTol', 1e-10);
        end

        function perturbed_trajectory_breaks_continuity(tc)
            nmpc = build_basic_nmpc();
            x_init = [0.5 -0.3];
            u_init = [0.1 0.2];

            w0 = nmpc.guess_from_initial(x_init, u_init);
            % Perturb one interior state in physical units and repack.
            [x, u] = nmpc.unpack_phys(w0);
            x(3, 1) = x(3, 1) + 1.0;
            ws = nmpc.pack_phys(x, u, []);
            [~, ceq] = nmpc.confun(ws, x_init);

            tc.verifyGreaterThan(norm(ceq, Inf), 1e-3);
        end

        function continuity_unaffected_by_x_scale(tc)
            % Continuity should hold under any non-trivial scaling.
            for sx = {[1 1], [10 0.1]}
                nmpc = build_basic_nmpc(x_scale=sx{1});
                x_init = [0.5 -0.3];
                u_init = [0.1 0.2];
                w0 = nmpc.guess_from_initial(x_init, u_init);
                [~, ceq] = nmpc.confun(w0, x_init);
                tc.verifyEqual(norm(ceq, Inf), 0, 'AbsTol', 1e-10);
            end
        end
    end
end
