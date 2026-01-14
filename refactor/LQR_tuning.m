clear all; close all; clc
rng(1)
%%
get_par

model = @(x, u)  dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

%% Linearization conditions
V = 1;
X = 10;

% Check for controllability
[xss, uss] = find_ss(V, X, par, model, ode_opt);
[A, B] = linearize(xss, uss, model);
controllability = ctrb(A,B);
unco = length(A) - rank(controllability);

%model = @(x, u) A*(x(:) - xss(:)) + B*(u(:) - uss(:));

%% Construct the LQR

% Targets
Xsp = 20;
Vsp = 1;
[xsp, usp] = find_ss(Vsp, Xsp, par, model, ode_opt);

% Sampling
Ts = 1/60;

% First try
Q = diag([1 1 1]);
R = diag([1 1 1]);
% Volume tightly
Q = diag([1 1 1]);
R = diag([1 1 1e-2]);
% % Volume and sugar tightly
% % Volume and sugar tightly
% Q = diag([100 1 1e3]);
% R = diag([10 1 1e-3]);



% Initial conditions
x0 = xss;
x0(2) = 2;
x0(3) = 1;

Tf = 20;

%%
var0 = [0    0   -1    -1   -1]; % Q(1,1) = 1
J0 = cost_LQR(var0, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);
opt_options = optimoptions('fminunc', 'Display','iter-detailed','UseParallel',true);
[optvar, J] = fminunc(@(var) cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt), var0, opt_options);

%%
K = build_LQR(optvar, A, B, Ts);

[Yode, Tode, U] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

t_interp = 0 : Ts : Tode(end);
X_interp = interp1(Tode, Yode(:, 2), t_interp);
Jx = trapz(t_interp, Xsp - X_interp);


colors = good_colors(3);
figure(2); clf
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

for i = 1:3
    nexttile
    plot(Tode, Yode(:,i), '-', 'LineWidth', 3, 'Color', colors(i,:)); hold on
    yline(xsp(i), '--', 'LineWidth', 3, 'Color', 'k')
    grid on; box on
    ylabel(sprintf('x_%d', i))
    if i == 1
        title('States')
        legend({'state','setpoint'}, 'Location','best', 'Interpreter','latex')
    end
    if i == 3
        xlabel('Time (h)')
    end
end

set_font_size()
%%
figure;
stairs(Tode, U)
function [Y, T, U] = LQR_simulation(system, Ts, Tf, y0, K, yss, uss, ode_opt)
% LQR_simulation for incremental LQR (delta-u LQR)
%
% Assumes K is designed for the augmented state:
%   z = [x_tilde; u_tilde_prev]
% where
%   x_tilde = x - yss
%   u_tilde_prev = (u_prev - uss)
% and the control law is
%   delta_u_tilde = -K * z
% so that
%   u = u_prev + delta_u_tilde
%
% Inputs:
%   system : function handle f(t,x,u) returning xdot (or x_{k+1} model)
%   Ts     : sample time (hours)
%   Tf     : final time (hours)
%   y0     : initial state (absolute)
%   K      : incremental LQR gain (size m x (n+m))
%   yss    : setpoint state (absolute)
%   uss    : setpoint input (absolute)
%   ode_opt: odeset options
%
% Outputs:
%   Y : state trajectory (absolute)
%   T : time vector
%   U : input trajectory (absolute), piecewise-constant over [T(k-1),T(k)]

    % Steps
    num_sim = ceil(Tf/Ts);

    % Dimensions
    n = length(yss);
    m = length(uss);

    % Preallocate
    Y = zeros(num_sim, n);
    U = zeros(num_sim, m);
    T = zeros(num_sim, 1);

    % Initialise
    Y(1,:) = y0(:).';
    T(1)   = 0;

    y_current = y0(:);
    t_current = 0;

    % Initialise previous input at setpoint (piecewise-constant hold)
    u_prev = uss(:);

    % (Optional) store also the initial input
    U(1,:) = u_prev(:).';

    for k = 2:num_sim
        % Deviations
        x_tilde      = y_current - yss(:);
        u_tilde_prev = u_prev    - uss(:);

        % Augmented deviation state
        z = [x_tilde; u_tilde_prev];

        % Incremental control in deviation coordinates
        du_tilde = -K * z;

        % Absolute input update
        uk = u_prev + du_tilde;

        % Constraints (keep exactly as you had them)
        uk = max(uk, zeros(m,1));
        uk = min(uk, 0.4*ones(m,1));
        if uk(1) > 0.4
            uk(1) = 0.4;
        end

        % ODE over the interval with piecewise-constant uk
        ode_current = @(t, y) system(t, y, uk);
        [t_span, y_span] = ode45(ode_current, [t_current, t_current + Ts], y_current, ode_opt);

        % Update state/time
        y_current = y_span(end,:).';
        y_current(y_current < 0) = 0;  % keep your non-negativity safeguard
        t_current = t_span(end);

        % Store
        Y(k,:) = y_current.';
        U(k,:) = uk.';
        T(k)   = t_current;

        % Update previous input for next increment
        u_prev = uk;
    end
end

function [xss, uss] = find_ss(V, X, par, model, ode_opt)
    % Starvation:
    % Solve for low flow rate steady-state
    Fin = X * par.Y_XSinv * par.kd / par.Sin; % approximate solution

    var = [Fin 0]; % solve ss for Fin and S
    var_opt = fsolve(@(var) model([V, X, var(2)], [var(1), 0, var(1)]), var);
    Fin = var_opt(1);
    S = var_opt(2);

    uss = [Fin 0 Fin];
    u = uss;
    TSPAN = [0 10];
    [t, x] = ode45(@(t,x) model(x, u), TSPAN, [V, X, S], ode_opt); % check
    xss = x(end, :);
end

function J = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt)
    K = build_LQR(var, A, B, Ts);

    [Yode, ~, Uode] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

    r = Yode - xsp;
    sr = sum(r.^2);
    ssr = sum(sr * [10; 1; 1]);
    reg = sum(sumsqr(diff(Uode)));
    J = ssr + reg;
end

function K = build_LQR(var, A, B, Ts)
    var = 10.^var;
    q = var(1:2);
    r = var(3:end);

    Q = diag([1 q]);
    R = diag(r);

    % LQR
    n = size(A,1);
    m = size(B,2);

    % Discretise
    sysc = ss(A,B,eye(n),zeros(n,m));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;
    
    % incremental form
    Ai = [Ad, Bd;
          zeros(m,n), eye(m)];
    Bi = [Bd;
          eye(m)];
    
    
    Qi = blkdiag(Q, zeros(m));  % Q only on the natural states
    Ri = R;                     % R acts on delta u explicitly
    
    K = dlqr(Ai, Bi, Qi, Ri);


    % The following yields an interesting structure, but cant be trusted
    % Ai = [A, B;
    %       zeros(m,n), zeros(m)];
    % Bi = [B;
    %       eye(m)];
    % Qi = blkdiag(Q, zeros(m));   % penalise x only
    % Ri = R*Ts^2;                 % The input is dudt and scalling matters
    % [K2,S,e] = lqrd(Ai,Bi,Qi,Ri,Ts);

end