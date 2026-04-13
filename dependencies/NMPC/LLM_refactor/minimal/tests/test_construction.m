classdef test_construction < matlab.unittest.TestCase
    % Construction, validation and default fill behaviour.

    methods (Test)
        %% Happy path
        function builds_with_defaults(tc)
            nmpc = build_basic_nmpc();
            tc.verifyEqual(nmpc.ny, nmpc.nx, 'ny defaults to nx');
            tc.verifyEqual(nmpc.x_scale, ones(1, nmpc.nx));
            tc.verifyEqual(nmpc.u_scale, ones(1, nmpc.nu));
            tc.verifyTrue(isempty(nmpc.P));
            tc.verifyTrue(isempty(nmpc.S));
            tc.verifyEqual(nmpc.rho_L1, 0);
        end

        function name_value_syntax_works(tc)
            % This exercises the modern foo(Name=Value) form.
            nmpc = build_basic_nmpc(P=eye(2));
            tc.verifyEqual(size(nmpc.P), [2 2]);
        end

        %% Required-field validation (cross-checked in init)
        function missing_required_errors(tc)
            % Omit `f`. Name-value arguments are always optional in MATLAB,
            % so the cross-check in init() is what raises here.
            tc.verifyError(@() NMPC(nx=2, nu=2, Ts=0.1, p=5, m=2, ...
                                    y_sp=[0 0], u_sp=[0 0], ...
                                    Q=eye(2), R=eye(2), ...
                                    Ymin=[-1 -1], Ymax=[1 1], ...
                                    umin=[-1 -1], umax=[1 1]), ...
                'NMPC:missing');
        end

        function unknown_field_errors(tc)
            tc.verifyError(@() build_basic_nmpc(nonsense_field=1), ...
                'MATLAB:TooManyInputs');
        end

        %% Cross-field validation (from init)
        function m_greater_than_p_errors(tc)
            tc.verifyError(@() build_basic_nmpc(p=3, m=5), 'NMPC:horizon');
        end

        function bad_Q_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Q=eye(3)), 'NMPC:Qsize');
        end

        function bad_R_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(R=eye(3)), 'NMPC:Rsize');
        end

        function bad_Ybounds_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Ymin=[0 0 0]), 'NMPC:Ybounds');
        end

        function bad_x_scale_errors(tc)
            tc.verifyError(@() build_basic_nmpc(x_scale=[1 -1]), 'NMPC:x_scale');
        end

        function bad_dumax_errors(tc)
            tc.verifyError(@() build_basic_nmpc(dumax=[1 0]), 'NMPC:dumax');
        end

        %% Slack mask handling
        function soft_mask_populates_indices(tc)
            nmpc = build_basic_nmpc(soft_mask=[true false], rho_L1=10);
            tc.verifyEqual(nmpc.soft_idx, 1);
            tc.verifyEqual(nmpc.n_soft, 1);
            tc.verifyEqual(nmpc.len_s, (nmpc.p + 1) * 1);
        end

        function no_soft_mask_means_no_slacks(tc)
            nmpc = build_basic_nmpc();
            tc.verifyEqual(nmpc.n_soft, 0);
            tc.verifyEqual(nmpc.len_s, 0);
        end

        %% solve() input validation
        function solve_rejects_wrong_size_x_init(tc)
            nmpc = build_basic_nmpc();
            tc.verifyError(@() nmpc.solve([0 0 0], [0 0]), 'NMPC:x_init');
        end

        function solve_rejects_wrong_size_u_init(tc)
            nmpc = build_basic_nmpc();
            tc.verifyError(@() nmpc.solve([0 0], [0 0 0]), 'NMPC:u_init');
        end

        function solve_rejects_nonfinite(tc)
            nmpc = build_basic_nmpc();
            tc.verifyError(@() nmpc.solve([NaN 0], [0 0]), ...
                'MATLAB:validators:mustBeFinite');
        end
    end
end
