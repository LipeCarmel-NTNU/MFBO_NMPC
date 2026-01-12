clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% Model
get_par
model = @(x, u)  dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-10, 'AbsTol', 1e-10);

%% Linearization conditions
V = 1;
X = 10;

% Check for controllability
[xss, uss] = find_ss(V, X, par, model, ode_opt);
[A, B] = linearize(xss, uss, model);
controllability = ctrb(A,B);
unco = length(A) - rank(controllability);

%% Construct the LQR

% Targets
Xsp = 20;
Vsp = 1;
[xsp, usp] = find_ss(Vsp, Xsp, par, model, ode_opt);

% Sampling
Ts = 20/3600;

% First test
% Q = diag([1 1 1]);
% R = diag([1 1 1]);

% Volume tightly controlled
% Q = diag([1 1 1]);
% R = diag([1 1 1e-2]);

% % Volume and sugar tightly controlled
opt_tuning = [-0.3628    4.5735   -3.9818    1.9878   -2.7478];
Q = diag([1 10.^opt_tuning(1:2)]);
R = diag(10.^opt_tuning(3:end));

[K,S,e] = lqrd(A,B,Q,R,Ts);

% Initial conditions
x0 = [1 2 0];

Tf = 20;
[Y, T, U] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

%%
t_interp = 0 : Ts : T(end);
X_interp = interp1(T, Y(:, 2), t_interp);
Jx = trapz(t_interp, Xsp - X_interp);

%%
colors = good_colors(3);
figure(2); %clf
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

for i = 1:3
    nexttile
    plot(T, Y(:,i), '-', 'LineWidth', 3, 'Color', colors(i,:)); hold on
    yline(xsp(i), '--', 'LineWidth', 3, 'Color', 'k')
    grid on; box on
    ylabel(sprintf('x_%d', i))
    if i == 1
        title('States')
        legend({'States','Setpoints'}, 'Location','best', 'Interpreter','latex')
    end
    if i == 3
        xlabel('Time (h)')
    end
end

set_font_size()

function [Y, T, U] = LQR_simulation(system, Ts, Tf, y0, K, yss, uss, ode_opt)
    % Number of steps
    num_sim = ceil(Tf/Ts);
    
    % Initialize output arrays
    ny = length(yss);        % Number of state variables
    nu = length(uss);        % Number of input variables
    Y = zeros(num_sim, ny);
    U = zeros(num_sim, nu); 
    T = zeros(num_sim, 1);  % Time vector
    
    % Set initial conditions
    Y(1, :) = y0;
    T(1) = 0;
    
    % Current state and time
    y_current = y0;
    t_current = 0;
    
    % Loop through each time step
    for k = 2 : num_sim
        % Measure
        ym = y_current(:) - yss(:);
        % ym(3) = max(ym(3) - 2, 0);
        %ym(1) = ym(1) + normrnd(0, 0.1);

        % Define the piecewise constant value of u
        uk = uss(:) - K*ym;
        uk = max(uk, zeros(nu, 1));
        uk = min(uk, 0.4*ones(nu, 1));
        if uk(1)>0.4
            uk(1)=0.4;
        end
        %disp(uk)

        % Define the ODE function for the current interval
        ode_current = @(t, y) system(t, y, uk);
        
        % Solve the ODE from t_current to t_current + Ts
        [t_span, y_span] = ode45(ode_current, [t_current, t_current + Ts], y_current, ode_opt);
        
        % Update current state and time
        y_current = y_span(end, :);
        y_current(y_current<0) = 0;
        t_current = t_span(end);
        
        % Store results
        Y(k, :) = y_current;
        U(k, :) = uk';
        T(k) = t_current;
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

