classdef test_construction < matlab.unittest.TestCase
    % Construction, validation and default fill behaviour.

    methods (Test)
        %% Happy path
        function builds_with_defaults(tc)
            nmpc = build_basic_nmpc();
            tc.verifyEqual(nmpc.x_scale, ones(1, nmpc.nx));
            tc.verifyEqual(nmpc.u_scale, ones(1, nmpc.nu));
            tc.verifyTrue(isempty(nmpc.P));
            tc.verifyTrue(isempty(nmpc.R_du));
            tc.verifyEqual(nmpc.rho_L1, 0);
        end

        function legacy_R_alias_sets_R_u(tc)
            nmpc = build_basic_nmpc(R_u=[], R=eye(2));
            tc.verifyEqual(nmpc.R_u, eye(2));
        end

        function conflicting_R_u_and_R_errors(tc)
            tc.verifyError(@() build_basic_nmpc(R_u=eye(2), R=2*eye(2)), ...
                'NMPC:R_u_alias');
        end

        function legacy_S_alias_sets_R_du(tc)
            nmpc = build_basic_nmpc(S=eye(2));
            tc.verifyEqual(nmpc.R_du, eye(2));
        end

        function conflicting_R_du_and_S_errors(tc)
            tc.verifyError(@() build_basic_nmpc(R_du=eye(2), S=2*eye(2)), ...
                'NMPC:R_du_alias');
        end

        function legacy_y_sp_alias_sets_x_sp(tc)
            nmpc = build_basic_nmpc(x_sp=[], y_sp=[3 4]);
            tc.verifyEqual(nmpc.x_sp, [3 4]);
        end

        function conflicting_x_sp_and_y_sp_errors(tc)
            tc.verifyError(@() build_basic_nmpc(x_sp=[1 2], y_sp=[3 4]), ...
                'NMPC:x_sp_alias');
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
                                    x_sp=[0 0], u_sp=[0 0], ...
                                    Q=eye(2), R_u=eye(2), ...
                                    Xmin=[-1 -1], Xmax=[1 1], ...
                                    umin=[-1 -1], umax=[1 1]), ...
                'NMPC:missing');
        end

        function unknown_field_errors(tc)
            tc.verifyError(@() build_basic_nmpc(nonsense_field=1), ...
                'MATLAB:TooManyInputs');
        end

        function output_map_is_not_configurable(tc)
            tc.verifyError(@() build_basic_nmpc(h_y=@(x) x), ...
                'MATLAB:TooManyInputs');
        end

        function integrator_is_not_configurable(tc)
            tc.verifyError(@() build_basic_nmpc(integrator=@(f, x, u) x), ...
                'MATLAB:TooManyInputs');
        end

        function terminal_setpoint_is_not_configurable(tc)
            tc.verifyError(@() build_basic_nmpc(x_sp_terminal=[0 0]), ...
                'MATLAB:TooManyInputs');
        end

        %% Cross-field validation (from init)
        function m_greater_than_p_errors(tc)
            tc.verifyError(@() build_basic_nmpc(p=3, m=5), 'NMPC:horizon');
        end

        function bad_Q_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Q=eye(3)), 'NMPC:Qsize');
        end

        function bad_R_u_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(R_u=eye(3)), 'NMPC:R_usize');
        end

        function bad_R_du_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(R_du=eye(3)), 'NMPC:R_dusize');
        end

        function legacy_Y_bounds_aliases_set_X_bounds(tc)
            nmpc = build_basic_nmpc(Xmin=[], Xmax=[], Ymin=[-3 -4], Ymax=[3 4]);
            tc.verifyEqual(nmpc.Xmin, [-3 -4]);
            tc.verifyEqual(nmpc.Xmax, [3 4]);
        end

        function conflicting_Xmin_and_Ymin_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Xmin=[-1 -2], Ymin=[-3 -4]), ...
                'NMPC:Xmin_alias');
        end

        function conflicting_Xmax_and_Ymax_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Xmax=[1 2], Ymax=[3 4]), ...
                'NMPC:Xmax_alias');
        end

        function bad_Xbounds_errors(tc)
            tc.verifyError(@() build_basic_nmpc(Xmin=[0 0 0]), 'NMPC:Xbounds');
        end

        function bad_x_scale_errors(tc)
            tc.verifyError(@() build_basic_nmpc(x_scale=[1 -1]), 'NMPC:x_scale');
        end

        function bad_dumax_errors(tc)
            tc.verifyError(@() build_basic_nmpc(dumax=[1 0]), 'NMPC:dumax');
        end

        function bad_rho_size_errors(tc)
            tc.verifyError(@() build_basic_nmpc(soft_mask=[true false], rho_L2=[1 2 3]), ...
                'NMPC:rho_size');
        end

        function rho_vector_must_match_nx_not_n_soft(tc)
            tc.verifyError(@() build_three_state_nmpc(rho_L2=[1 1]), ...
                'NMPC:rho_size');
        end

        function rho_vector_with_length_nx_is_valid(tc)
            nmpc = build_three_state_nmpc(rho_L2=[1 1 0]);
            tc.verifyEqual(nmpc.rho_L2, [1 1 0]);
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

function nmpc = build_three_state_nmpc(varargin)
    defaults = struct( ...
        'xdot',      @(x, u) [-x(1) + u(1), -2*x(2) + u(2), -0.5*x(3)], ...
        'nx',        3, ...
        'nu',        2, ...
        'Ts',        0.1, ...
        'p',         4, ...
        'm',         2, ...
        'x_sp',      [0 0 0], ...
        'u_sp',      [0 0], ...
        'Q',         eye(3), ...
        'R_u',       eye(2), ...
        'Xmin',      [-1 -1 -1], ...
        'Xmax',      [1 1 1], ...
        'umin',      [-1 -1], ...
        'umax',      [1 1], ...
        'soft_mask', [true true false]);

    overrides = struct(varargin{:});
    fn = fieldnames(overrides);
    for k = 1:numel(fn)
        defaults.(fn{k}) = overrides.(fn{k});
    end

    args = namedargs2cell(defaults);
    nmpc = NMPC(args{:});
end
