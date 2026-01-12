% Casadi implementation of orthogonal collocation adapted to implement flexible step length
%
%     Original code:
%     https://github.com/casadi/casadi/blob/main/docs/examples/matlab/direct_collocation.m
%     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
%
%     For the modified collocation see page 297 of:
%     Biegler, L.T. (2010). Nonlinear Programming: Concepts,
%     Algorithms, and Applications to Chemical Processes.
%     MOS-SIAM Series on Optimization. Society for Indus-
%     trial and Applied Mathematics : Mathematical Pro-
%     gramming Society, Philadelphia.


% An implementation of direct collocation
% Joel Andersson, 2016
clearvars; close all; clc
import casadi.*

u_ref = 0.0; % Fails for low D

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time horizon
T = 10;
x0_init = [10; 3];
N = 50; % number of discretization intervals

model = @ monod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
u = SX.sym('u');

% Model equations
xdot = model(0, [x1; x2], u);

% Objective term
L = 0*(x1^2 + x2^2) + (u-u_ref)^2;

% Continuous time dynamics
f = Function('f', {x, u}, {xdot, L});

% Control discretization
h = T/N;

%--------------------------
% Predefine bounds
%--------------------------
% States: x = [x1; x2]
lb_x = [0; 0];   % lower bound for (x1, x2)
ub_x = [ inf;  inf];    % upper bound for (x1, x2)
x0_guess = x0_init;      % initial guess for (x1, x2)

% Control: u
lb_u = 0;
ub_u = 0.5;
u0_guess = u_ref;

% Initial condition

%--------------------------
degree = 3; method = 'radau';
[C, D, B] = collocation_polynomial(degree, method);

%% --------------------------
% Pre-initialize symbolic variables as cells
%--------------------------
% State variables at interval boundaries (N+1 total: 0, 1, ..., N)
X = cell(N+1, 1);
for k = 0:N
    X{k+1} = MX.sym(['X_' num2str(k)], 2);
end

% Control variables (N total: 0, 1, ..., N-1)
U = cell(N, 1);
for k = 0:N-1
    U{k+1} = MX.sym(['U_' num2str(k)]);
end

% Collocation state variables (N intervals × degree collocation points)
X_coll = cell(N, degree);
for k = 0:N-1
    for j = 1:degree
        X_coll{k+1, j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2);
    end
end

% Flexible h
hk = cell(N);
for k = 1:N
    hk{k} = MX.sym(['h_' num2str(k)]);
end


% --------------------------
% Build constraints and objective (with step regularisation)
% --------------------------
J = 0;
g = {};
lbg = [];
ubg = [];

h_nom = T/N;        % nominal step size
gamma = 0.9;       % step variability bound
rho_h  = 1e-1;      % small regularisation weight on (hk - h_nom)^2

% Formulate the NLP
for k = 0:N-1
    k_idx = k + 1;

    % Current step length for this interval
    hk_k = hk{k_idx};

    % Collocation states at the end of the interval
    Xk_end = D(1) * X{k_idx};

    for j = 1:degree
        % Expression for the state derivative at the collocation point
        xp = C(1, j+1) * X{k_idx};
        for r = 1:degree
            xp = xp + C(r+1, j+1) * X_coll{k_idx, r};
        end

        % Append collocation equations
        [fj, qj] = f(X_coll{k_idx, j}, U{k_idx});

        % Collocation (defect) constraint: xp = hk * f
        g   = {g{:}, hk_k * fj - xp};
        lbg = [lbg; 0; 0];
        ubg = [ubg; 0; 0];

        % Add contribution to the end state
        Xk_end = Xk_end + D(j+1) * X_coll{k_idx, j};

        % Add contribution to quadrature function
        J = J + B(j+1) * qj * hk{k_idx};
    end

    % Add equality constraint for continuity
    g = {g{:}, Xk_end - X{k_idx + 1}};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];

    % --- Step regularisation (keeps hk close to nominal) ---
    J = J + rho_h * (hk_k - h_nom)^2;

end

% --------------------------
% Compose decision variables and bounds
% --------------------------
w = {};
w0 = [];
lbw = [];
ubw = [];

% Add state variables at interval boundaries
for k = 0:N
    k_idx = k + 1;
    w = {w{:}, X{k_idx}};

    if k == 0
        lbw = [lbw; x0_init];
        ubw = [ubw; x0_init];
        w0  = [w0; x0_init];
    else
        lbw = [lbw; lb_x];
        ubw = [ubw; ub_x];
        w0  = [w0; x0_guess];
    end
end

% Add control variables
for k = 0:N-1
    k_idx = k + 1;
    w = {w{:}, U{k_idx}};
    lbw = [lbw; lb_u];
    ubw = [ubw; ub_u];
    w0  = [w0; u0_guess];
end

% Add collocation state variables
for k = 0:N-1
    k_idx = k + 1;
    for j = 1:degree
        w = {w{:}, X_coll{k_idx, j}};
        lbw = [lbw; lb_x];
        ubw = [ubw; ub_x];
        w0  = [w0; x0_guess];
    end
end

% Add step size variables + bounds + nominal initial guess
sum_h = 0;
for k = 1:N
    w = {w{:}, hk{k}};
    lbw = [lbw; h_nom*(1 - gamma)];
    ubw = [ubw; h_nom*(1 + gamma)];
    w0  = [w0; h_nom];

    sum_h = sum_h + hk{k};
end

% Enforce total horizon exactly: sum(hk) = T
g   = {g{:}, sum_h};
lbg = [lbg; T];
ubg = [ubg; T];


%%
clc
% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

% Solve the NLP

solver = nlpsol('solver', 'ipopt', prob);
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);


% solver = nlpsol('solver', 'sqpmethod', prob);
% sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
%             'lbg', lbg, 'ubg', ubg);
% w_opt = full(sol.x);

%% --------------------------
% Extract solution
%--------------------------
% Calculate indices for solution extraction
n_states = N + 1;  % Number of state variables at boundaries
n_controls = N;    % Number of control variables
n_coll_states = N * degree * 2;  % Number of collocation state variables (2 states each)

u_idx = n_states*2;
coll_idx = u_idx + n_controls;
h_idx = coll_idx + n_coll_states;

% Extract states at interval boundaries
x_states = reshape(w_opt(1:n_states*2), 2, n_states)';  % Each state is 2D
x1_opt = x_states(:, 1);
x2_opt = x_states(:, 2);

% Extract controls
u_opt = w_opt(u_idx + 1 : n_states*2 + n_controls);

t_opt = [0; cumsum(w_opt(h_idx + 1 : h_idx + N))];

%% Reference solution
clf;
hold on
opts = odeset('NonNegative', 2, 'RelTol', 1e-10, 'AbsTol', 1e-10);
[t_ode, x_ode] = ode45(@(t,y) model(t,y, u0_guess), [0 T], x0_init, opts);
plot(t_ode, x_ode, 'k--')

%% Plot the solution
x_coll = [x1_opt, x2_opt];

%tgrid = linspace(0, T, N+1);
tgrid = t_opt;
plot(tgrid, x1_opt, 'bo')
plot(tgrid, x2_opt, 'ro')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
%legend('x1','x2','u')

%% Check error
x_ode_ref = interp1(t_ode, x_ode, tgrid);
e = x_coll - x_ode_ref;

SSE = sum(sumsqr(e)); % 1.2546e-09
MaxAE = max(max(abs(e))); % 1.0793e-05

disp('SSE and Max abs error:')
disp(SSE)
disp(MaxAE)
% if SSE > 2e-9 || MaxAE > 2e-5
%     error('Check the solution')
% end

function dxdt = simple_model(t,y,u)
    x1 = y(1);
    x2 = y(2);
    dxdt = [- x1; -x2*x1 - x2*0.1];
end

function dxdt = monod(t,x,u)
    par = struct();
    
    par.Sin                    = 70;
    par.KcS                    = 0.0321;
    par.tau_iS                 = 4.2100;
    par.Y_XS                   = 0.4204;
    par.Y_XSinv                = 2.3786;
    par.mu_max                 = 0.1945;
    par.Ks                     = 0.0070;
    par.Y_CO2X                 = 0.5430;
    par.kd                     = 0.0060;
    par.Y_CO2Xinv               = 0.3581;
    par.mu                     = 0.1945;
    % States
    X = x(1); S = x(2); D = u;
    
    Sin= par.Sin;
    
    % Parameters
    Y_XS   = par.Y_XS;
    mu_max = par.mu_max;    % h^-1
    Ks     = par.Ks;        % g L^-1 
    kd     = par.kd;
    
    % Define the rate
    mu = mu_max .*(S ./(Ks + S));
    
    % The inputs are Fin, Fout
    % Differential equations:
    dX   = -X.*(D) +  mu.*X - kd .*X; %Biomass
    dS   = (Sin-S).*(D) - mu .*X ./Y_XS; % Substrate
    
    % Output:
    dxdt = [dX; dS];
end
function [C, D, B] = collocation_polynomial(degree, method)
    import casadi.*
    
    % Get collocation points
    tau_root = [0 collocation_points(degree, method)];
    
    % Coefficients of the collocation equation
    C = zeros(degree+1,degree+1);
    
    % Coefficients of the continuity equation
    D = zeros(degree+1, 1);
    
    % Coefficients of the quadrature function
    B = zeros(degree+1, 1);
    
    % Construct polynomial basis
    for j=1:degree+1
      % Construct Lagrange polynomials to get the polynomial basis at the collocation point
      coeff = 1;
      for r=1:degree+1
        if r ~= j
          coeff = conv(coeff, [1, -tau_root(r)]);
          coeff = coeff / (tau_root(j)-tau_root(r));
        end
      end
      % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D(j) = polyval(coeff, 1.0);
    
      % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      pder = polyder(coeff);
      for r=1:degree+1
        C(j,r) = polyval(pder, tau_root(r));
      end
    
      % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      pint = polyint(coeff);
      B(j) = polyval(pint, 1.0);
    end
end

