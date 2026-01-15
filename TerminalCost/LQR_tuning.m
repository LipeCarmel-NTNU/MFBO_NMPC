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
Xsp = 15;
Vsp = 1;
[xsp, usp] = find_ss(Vsp, Xsp, par, model, ode_opt);

% Sampling
Ts = 1/60;

% Initial conditions
x0 = xss;
x0(1) = 1.1;
x0(2) = 5;
x0(3) = 1;

Tf = 20;

%% Tuning
var0 =  [0 0 0 0 log10(0.5)]; % Q(1,1) = 1
J0 = cost_LQR(var0, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);
opt_options = optimoptions('fminunc', 'Display','iter-detailed','UseParallel',true, 'MaxFunctionEvaluations', 2000);
[optvar, J] = fminunc(@(var) cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt), var0, opt_options);

%% Simulation
tuning_par = optvar;
[K, Qi, Ri] = build_LQR(tuning_par, A, B, Ts);

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

figure(1);
stairs(Tode, U)
figure(2)



function J = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt)
    K = build_LQR(var, A, B, Ts);

    [Yode, ~, Uode] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

    r = Yode - xsp;
    r = r.*[10 1 1];
    sr = sum(r.^2);
    ssr = sum(sr);                             % 1.7539e+03
    reg = 1e4*sum(sumsqr(diff(Uode)));         % 173.2546        1e.2*sum
    reg2 = 1e2*sumsqr(var);                    %  25.4240
    J = ssr + reg + reg2;
    disp(J)
    disp(reg2)
end
