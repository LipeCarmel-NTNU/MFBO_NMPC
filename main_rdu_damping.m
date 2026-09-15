%% MAIN_RDU_DAMPING  Damp MF_1 towards BO_1's input activity, keeping structure.
%
% The blend is a single scalar: Rdu = rho * Rdu_MF, rho >= 1. Scaling keeps
% the RELATIVE weighting between the three input channels that the BO found,
% and only changes how hard the controller is pushed overall. Q, Np and Nc
% stay at MF_1's values, so "the same tracking trajectory" means the same
% objective, not a re-tuned one.
%
% rho is chosen on the LINEARISED, STOCHASTIC steady state. Near the setpoint
% the setpoint is constant and the state is already there, so what moves the
% input is not tracking but measurement noise: nmpc_run_case feeds
% yk_meas = xk + noise straight into the solver, with no filter. With
% Ru = 0 (disabled in both arms) the incremental LQR is
%
%   K  = dlqr(Ai, Bi, blkdiag(Q, 0), Rdu),      z = [x - xss; u_prev - uss]
%   Du = -K z - Kx v,                           v ~ N(0, W), W = diag(sigma_y^2)
%   z+ = (Ai - Bi K) z - Bi Kx v
%
% so the stationary covariance Sigma solves one discrete Lyapunov equation and
%
%   Jx(Rdu)  = tr( Q Sx Sigma Sx' )                 steady-state tracking
%   Jdu(Rdu) = tr( K Sigma K' ) + tr( Kx W Kx' )    steady-state input motion
%
% The program solved here is
%
%   min over rho >= 1   Jdu(rho * Rdu_MF)
%   s.t.  Jx (rho * Rdu_MF) >= (1 - eps) * Jx(Rdu_MF)      retain the loop gain
%         Jdu(rho * Rdu_MF) >= Jdu(Rdu_BO)                 do not damp past BO_1
%         rho * Rdu_MF      <= 10^THETA_MAX                stay in the BO box
%
% Jx and Jdu are BOTH decreasing in rho, so the minimisation drives rho up and
% every constraint is an upper bound; the answer is the smallest of the three.
%
% The Jx constraint is a FLOOR, not a budget, and the direction is deliberate.
% Jx is the stationary state covariance driven by measurement noise, and that
% noise reaches the state only through K: with K = 0 the linearised plant sits
% at the setpoint forever and Jx is exactly zero. Jx therefore decreases
% monotonically to zero in rho - the do-nothing controller is its exact
% minimiser - so an UPPER budget on Jx would be satisfied by the argmin of the
% objective and could never bind. Holding Jx at (1 - eps) of MF_1's value
% instead asks the loop to keep 90% of the gain it had. That binds, and it is
% the knob.
%
% Name it accordingly: Jx is RETAINED LOOP GAIN, or closed-loop sensitivity to
% measurement noise. It is NOT disturbance rejection. base.plant = base.model
% and the only stochastic input in the campaign is measurement noise, so there
% is no process disturbance here to reject, and none is invented to make the
% constraint bite.
%
% Jx_tr and Jdu_tr are printed beside them as a zero-parameter cross-check:
% the deterministic closed-loop cost from the two initial conditions the
% campaign actually ran, under the diag([10 1 1])^2 weighting finalize_case
% already applies. They assume nothing that is not already in the experiment.
% At eps = 0.10 the two measures agree closely, which is the argument for
% steering by the cheap one.
%
% Nothing here is written into a campaign folder. The results land in
% results/rdu_damping/ with no results.csv, so parse_registry.py will not
% discover it as a case.

clear; close all; clc;

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% ------------------------------------------------------------------ knobs
EPS_X       = 0.90;     % keep at least (1 - EPS_X) of MF_1's loop gain Jx
SIM_HOURS   = 10;        % horizon of the verification run
RNG_SEED    = 1;        % fixes base.noise, shared by every simulated theta
THETA_MAX   = 3;        % upper bound of the r_du exponent in the BO box
% Every simulation is expensive, so only the blend is run. MF_1's trajectory is
% already on disk and re-simulating it here buys only the identical noise
% realisation; BO_1 costs roughly 30 minutes on its own (Np = 15, Nc = 6).
SIM_MF      = false;
SIM_BO      = false;
X0          = [1.0, 10, 0; 1.1, 25, 5];   % the campaign's initial conditions

MF_FILE = fullfile(current_dir, 'results', 'case2_v3',        'out_20260911_164617.mat');
BO_FILE = fullfile(current_dir, 'results', 'results_baseline', 'out_20260907_155726.mat');
OUT_DIR = fullfile(current_dir, 'results', 'rdu_damping');
DAMPING_FILENAME = 'rdu_damping_EPS_90.mat';
%% ------------------------------------------------------- the two parents
assert(isfile(MF_FILE), 'MF_1 trends file not found: %s', MF_FILE);
assert(isfile(BO_FILE), 'BO_1 trends file not found: %s', BO_FILE);
Smf = load(MF_FILE, 'out');  theta_mf = Smf.out.theta(:).';
Sbo = load(BO_FILE, 'out');  theta_bo = Sbo.out.theta(:).';

rng(RNG_SEED)
base = nmpc_base();                 % same sigma_y, Ts, tf and setpoint as the campaign
nx = base.nx;  nu = base.nu;

cfg_mf = decode_theta(theta_mf, nx, nu);
cfg_bo = decode_theta(theta_bo, nx, nu);

assert(max(abs(diag(cfg_mf.Ru))) == 0 && max(abs(diag(cfg_bo.Ru))) == 0, ...
    ['Ru is nonzero in one of the parents. The synthesis below assumes the ' ...
     'disabled Ru of case2, which makes the LQR cross term vanish.']);

%% ------------------------------------- linearised incremental plant model
[A, B]   = local_linearize(base);
[Ai, Bi] = incremental(A, B, base.dt);

W  = diag(base.sigma_y(:).^2);
Sx = [eye(nx), zeros(nx, nu)];

% Zero-parameter cross-check. Wt is the error weighting finalize_case applies
% when it forms the reported J_track; Z0 is the campaign's two initial
% conditions in incremental coordinates, and nmpc_run_case starts from
% uk = zeros(1, nu), so u_prev = 0 rather than uss.
Wt = diag([10 1 1]).^2;
Z0 = [X0 - base.xsp, zeros(size(X0, 1), nu) - base.usp];

Q      = cfg_mf.Q;          % the weights that define "MF_1's trajectory"
Rdu_mf = cfg_mf.Rdu;
Rdu_bo = cfg_bo.Rdu;

%% ---------------------------------------- steady-state costs of the parents
[Jx_mf, Jdu_mf] = ss_costs(Ai, Bi, Q,       Rdu_mf, W, Sx, nx);
[Jx_bo, Jdu_bo] = ss_costs(Ai, Bi, cfg_bo.Q, Rdu_bo, W, Sx, nx);

% Jdu_bo is BO_1's own input activity, under BO_1's own Q. That is the level
% the blend is asked not to go below. The same Rdu read under MF_1's Q is a
% different number and is printed only for reference.
[~, Jdu_bo_underQmf] = ss_costs(Ai, Bi, Q, Rdu_bo, W, Sx, nx);

[Jxt_mf, Jdut_mf] = tr_costs(Ai, Bi, Q,        Rdu_mf, Z0, Wt, Sx);
[Jxt_bo, Jdut_bo] = tr_costs(Ai, Bi, cfg_bo.Q, Rdu_bo, Z0, Wt, Sx);

%% --------------------------------------------------------- solve for rho
rho_box = min(10^THETA_MAX ./ diag(Rdu_mf));
assert(rho_box > 1, 'Rdu_MF already sits at the top of the box; no room to damp.');

fJx  = @(r) ss_costs(Ai, Bi, Q, r * Rdu_mf, W, Sx, nx);
fJdu = @(r) second_out(@() ss_costs(Ai, Bi, Q, r * Rdu_mf, W, Sx, nx));

rho_eps = largest_feasible(fJx,  @(v) v >= (1 - EPS_X) * Jx_mf,  1, rho_box);
rho_flr = largest_feasible(fJdu, @(v) v >= Jdu_bo,               1, rho_box);

if isnan(rho_flr)
    error(['Infeasible at rho = 1: MF_1 already moves the inputs less than ' ...
           'BO_1 in the linear steady state (Jdu_MF = %.6g < Jdu_BO = %.6g). ' ...
           'There is nothing to damp.'], Jdu_mf, Jdu_bo);
end
if isnan(rho_eps)
    error('Infeasible at rho = 1: Jx(1) < (1-eps) Jx(1). Check 0 <= EPS_X < 1.');
end

[rho, which] = min([rho_eps, rho_flr, rho_box]);
names  = ["Jx floor (1-eps)", "Jdu floor at BO_1", "r_du box edge"];
active = names(which);

Rdu_new   = rho * Rdu_mf;
theta_new = theta_mf;
theta_new(3 + nx + nu + (1:nu)) = log10(diag(Rdu_new)).';   % the r_du block
[Jx_new, Jdu_new]   = ss_costs(Ai, Bi, Q, Rdu_new, W, Sx, nx);
[Jxt_new, Jdut_new] = tr_costs(Ai, Bi, Q, Rdu_new, Z0, Wt, Sx);

%% ------------------------------------------------------------ the printout
lab = ["MF_1", "BLEND", "BO_1"];
fprintf('\n================ Rdu damping synthesis ================\n');
fprintf('steady state linearised at xss = [%.4f %.4f %.6g], uss = [%.6g %.6g %.6g]\n', ...
    base.xsp, base.usp);
fprintf('sigma_y = [%g %g %g],  Ts = %g h\n\n', base.sigma_y, base.dt);

fprintf('rho*              = %.6g      (active constraint: %s)\n', rho, active);
fprintf('  rho from Jx floor   = %.6g   (keep >= %.1f%% of MF_1''s loop gain)\n', ...
    rho_eps, 100 * (1 - EPS_X));
fprintf('  rho from Jdu floor  = %.6g   (do not damp past BO_1)\n', rho_flr);
fprintf('  rho from box edge   = %.6g\n', rho_box);
if rho_eps >= rho_box
    fprintf(['  NOTE: the Jx floor did not bind before the box edge, so EPS_X is\n' ...
             '        loose enough that the detuning budget is not what limits rho.\n']);
end

fprintf('\n%-22s%14s%14s%14s\n', 'tuning', lab(1), lab(2), lab(3));
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-22s%14d%14d%14d\n', 'Np', cfg_mf.p, cfg_mf.p, cfg_bo.p);
fprintf('%-22s%14d%14d%14d\n', 'Nc', cfg_mf.m, cfg_mf.m, cfg_bo.m);
print_row('theta_q',   log10(diag(cfg_mf.Q)),   log10(diag(Q)),        log10(diag(cfg_bo.Q)));
print_row('theta_rdu', log10(diag(Rdu_mf)),     log10(diag(Rdu_new)),  log10(diag(Rdu_bo)));
print_row('Rdu',       diag(Rdu_mf),            diag(Rdu_new),         diag(Rdu_bo));
print_row('Rdu / Q(3)', diag(Rdu_mf)/Q(3,3),    diag(Rdu_new)/Q(3,3),  diag(Rdu_bo)/cfg_bo.Q(3,3));

fprintf('\n%-22s%14s%14s%14s\n', 'steady-state cost', lab(1), lab(2), lab(3));
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-22s%14.6g%14.6g%14.6g\n', 'Jx',            Jx_mf,  Jx_new,  Jx_bo);
fprintf('%-22s%14.6g%14.6g%14.6g\n', 'Jdu',           Jdu_mf, Jdu_new, Jdu_bo);
fprintf('%-22s%14.4f%14.4f%14.4f\n', 'Jx  / Jx(MF_1)',  1, Jx_new/Jx_mf,   Jx_bo/Jx_mf);
fprintf('%-22s%14.5f%14.5f%14.5f\n', 'Jdu / Jdu(MF_1)', 1, Jdu_new/Jdu_mf, Jdu_bo/Jdu_mf);
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-22s%14.6g%14.6g%14.6g\n', 'Jx_tr  (check)',  Jxt_mf,  Jxt_new,  Jxt_bo);
fprintf('%-22s%14.6g%14.6g%14.6g\n', 'Jdu_tr (check)',  Jdut_mf, Jdut_new, Jdut_bo);
fprintf('%-22s%14.4f%14.4f%14.4f\n', 'Jx_tr  / MF_1',   1, Jxt_new/Jxt_mf,   Jxt_bo/Jxt_mf);
fprintf('%-22s%14.4f%14.4f%14.4f\n', 'Jdu_tr / MF_1',   1, Jdut_new/Jdut_mf, Jdut_bo/Jdut_mf);
fprintf('\n(Jdu of BO_1''s Rdu read under MF_1''s Q, for reference: %.6g)\n', Jdu_bo_underQmf);

fprintf('\ntheta_blend = [');  fprintf(' %.6g', theta_new);  fprintf(' ]\n');
fprintf('======================================================\n\n');

pause(1)

%% ------------------------------------------------- verification simulation
if ~isfolder(OUT_DIR); mkdir(OUT_DIR); end
z_sim = SIM_HOURS / base.tf;
assert(z_sim > 0 && z_sim <= 1, 'SIM_HOURS must lie in (0, base.tf].');

runs = struct('label', {}, 'theta', {}, 'out', {});
runs(end+1) = struct('label', "BLEND", 'theta', set_z(theta_new, z_sim), 'out', []);
if SIM_MF; runs(end+1) = struct('label', "MF_1", 'theta', set_z(theta_mf, z_sim), 'out', []); end
if SIM_BO; runs(end+1) = struct('label', "BO_1", 'theta', set_z(theta_bo, z_sim), 'out', []); end

for k = 1:numel(runs)
    fprintf('simulating %s over %g h ...\n', runs(k).label, SIM_HOURS);
    t0 = tic;
    runs(k).out = simulate_nmpc(base, runs(k).theta, ...
        horizon = "fidelity", extrapolate = false, verbosity = "control", ...
        x0 = X0, run_id = "rdu_damping_" + runs(k).label, log_path = "");
    fprintf('   done in %.1f s wall, %.1f s solver\n', toc(t0), runs(k).out.runtime_s);
end

%% ------------------------------------------------------- measured results
stateNames = ["V", "X", "S"];
cols = [runs.label];
fprintf('\n============ measured over %g h, same noise realisation ============\n', SIM_HOURS);
fprintf('%-16s', 'metric'); fprintf('%14s', cols); fprintf('\n');
fprintf('%s\n', repmat('-', 1, 16 + 14*numel(cols)));
fprintf('%-16s', 'J_track'); fprintf('%14.2f', arrayfun(@(r) r.out.SSE,  runs)); fprintf('\n');
fprintf('%-16s', 'J_TV');    fprintf('%14.5f', arrayfun(@(r) r.out.SSdU, runs)); fprintf('\n');
fprintf('%-16s', 't_nmpc (s)'); fprintf('%14.1f', arrayfun(@(r) r.out.runtime_s, runs)); fprintf('\n');
fprintf('%s\n', repmat('-', 1, 16 + 14*numel(cols)));
for s = 1:numel(runs(1).out.case)
    for j = 1:nx
        fprintf('%-16s', sprintf('IAE_%s_%d', stateNames(j), s));
        fprintf('%14.5g', arrayfun(@(r) iae(r.out.case(s), j, base.dt), runs)); fprintf('\n');
    end
    for j = 1:nu
        fprintf('%-16s', sprintf('TV_u%d_%d', j, s));
        fprintf('%14.5g', arrayfun(@(r) tv(r.out.case(s), j), runs)); fprintf('\n');
    end
    fprintf('%s\n', repmat('-', 1, 16 + 14*numel(cols)));
end

synth = struct('rho', rho, 'active_constraint', active, 'EPS_X', EPS_X, ...
    'theta_mf', theta_mf, 'theta_bo', theta_bo, 'theta_blend', theta_new, ...
    'Jx', [Jx_mf Jx_new Jx_bo], 'Jdu', [Jdu_mf Jdu_new Jdu_bo], ...
    'Jx_tr', [Jxt_mf Jxt_new Jxt_bo], 'Jdu_tr', [Jdut_mf Jdut_new Jdut_bo], ...
    'Ai', Ai, 'Bi', Bi, 'W', W, 'sim_hours', SIM_HOURS, 'rng_seed', RNG_SEED);
save(fullfile(OUT_DIR, DAMPING_FILENAME), 'synth', 'runs');
fprintf('\nwrote %s\n', fullfile(OUT_DIR, DAMPING_FILENAME));


%% ===================== local functions =====================

function [Jx, Jdu, K, Sigma] = ss_costs(Ai, Bi, Q, Rdu, W, Sx, nx)
%SS_COSTS Stationary tracking and input-motion cost of the incremental LQR.
%   Ru is zero, so the dlqr cross term N = [0; Ru] vanishes and R = Rdu.
    Qz = blkdiag(Q, zeros(size(Rdu)));
    K  = dlqr(Ai, Bi, Qz, Rdu);
    Kx = K(:, 1:nx);
    Acl = Ai - Bi * K;
    % Measurement noise enters the loop only through the controller, because
    % the solver is handed yk_meas = xk + v with no filter in between.
    Sigma = dlyap(Acl, Bi * Kx * W * Kx.' * Bi.');
    Jx  = trace(Q * Sx * Sigma * Sx.');
    % v_k is independent of z_k, so the two terms add.
    Jdu = trace(K * Sigma * K.') + trace(Kx * W * Kx.');
end

function [Jx, Jdu] = tr_costs(Ai, Bi, Q, Rdu, Z0, Wt, Sx)
%TR_COSTS Deterministic closed-loop cost from the campaign's initial states.
%   The zero-parameter cross-check on the steady-state pair: no noise model,
%   no invented disturbance, just the two initial conditions that were run and
%   the weighting finalize_case already uses.
    Qz  = blkdiag(Q, zeros(size(Rdu)));
    K   = dlqr(Ai, Bi, Qz, Rdu);
    Acl = Ai - Bi * K;
    Px  = dlyap(Acl.', Sx.' * Wt * Sx);
    Pdu = dlyap(Acl.', K.' * K);
    Jx  = 0;
    Jdu = 0;
    for i = 1:size(Z0, 1)
        z   = Z0(i, :).';
        Jx  = Jx  + z.' * Px  * z;
        Jdu = Jdu + z.' * Pdu * z;
    end
end

function rho = largest_feasible(fun, ok, lo, hi)
%LARGEST_FEASIBLE Biggest rho in [lo, hi] with ok(fun(rho)) true.
%   Assumes ok(fun(.)) is true on an interval starting at lo. Bisects on
%   log10(rho). NaN when even lo fails.
    if ~ok(fun(lo)); rho = NaN; return; end
    if  ok(fun(hi)); rho = hi;  return; end
    a = log10(lo); b = log10(hi);
    for i = 1:100
        c = 0.5 * (a + b);
        if ok(fun(10^c)); a = c; else; b = c; end
    end
    rho = 10^a;
end

function v = second_out(fh)
%SECOND_OUT The second output of a zero-argument handle, as a value.
    [~, v] = fh();
end

function [A, B] = local_linearize(base)
%LOCAL_LINEARIZE Jacobians at the setpoint, symbolic if available.
    try
        [A, B] = linearize(base.xsp, base.usp, base.model);
    catch ME
        warning('linearize failed (%s); falling back to central differences.', ME.message);
        A = fd_jac(@(x) base.model(x, base.usp), base.xsp);
        B = fd_jac(@(u) base.model(base.xsp, u), base.usp);
    end
end

function J = fd_jac(fh, z0)
%FD_JAC Central-difference Jacobian of fh at z0.
    z0 = z0(:).';
    f0 = fh(z0);
    n  = numel(z0);
    J  = zeros(numel(f0), n);
    for i = 1:n
        h = 1e-7 * max(1, abs(z0(i)));
        zp = z0; zp(i) = zp(i) + h;
        zm = z0; zm(i) = zm(i) - h;
        J(:, i) = (fh(zp) - fh(zm)) / (2 * h);
    end
end

function theta = set_z(theta, z)
%SET_Z Replace the fidelity component, leaving every weight alone.
    theta(1) = z;
end

function val = iae(c, j, dt_h)
%IAE Integral absolute error of state j, rectangular rule, true state.
    n = double(c.i_last);
    val = sum(abs(c.Y(1:n, j) - c.Ysp(1:n, j))) * dt_h;
end

function val = tv(c, j)
%TV Total variation of input j.
    n = double(c.i_last);
    val = sum(abs(diff(c.U(1:n, j))));
end

function print_row(name, a, b, c)
%PRINT_ROW One tuning row, three columns, one line per input or state.
    a = a(:); b = b(:); c = c(:);
    for i = 1:numel(a)
        if i == 1; tag = name; else; tag = ''; end
        fprintf('%-22s%14.6g%14.6g%14.6g\n', sprintf('%s(%d)', tag, i), a(i), b(i), c(i));
    end
end
