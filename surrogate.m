close all; clc
load('data.mat')

% -----------------------------
% Build cumulative signals
% -----------------------------
N = length(out.case(1).SSE);
f = (1:N).' / N;                 % p in (0,1]
x = 2*f - 1;                     % map to [-1,1]

time = cumsum(out.case(1).RUNTIME) + cumsum(out.case(2).RUNTIME);

SSdU = [0; cumsum(out.case(1).SSdU)] + [0; cumsum(out.case(2).SSdU)];
SSdU = SSdU / SSdU(end);

SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
SSE  = SSE / SSE(end);

% -----------------------------
% Fit with fminunc
% -----------------------------
opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','off', ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);

% Cheb5: coefficients c = [c0..c5]
c5_0 = zeros(6,1);

c_SSdU = fminunc(@(c) obj_Cheb5(c, x,  SSdU), c5_0, opts);
c_SSE  = fminunc(@(c) obj_Cheb5(c, x,  SSE ), c5_0, opts);

SSdU_fit = Cheb5(x, c_SSdU);
SSE_fit  = Cheb5(x, c_SSE);


c2_0 = zeros(4,1);

c_time = fminunc(@(c) obj_Cheb2_constrained(c, x, time), c2_0, opts);
time_fit = Cheb_time(x, c_time);

% -----------------------------
% Plots
% -----------------------------
figure;
plot(f, SSdU, 'LineWidth', 1.5); hold on
plot(f, SSdU_fit, '--', 'LineWidth', 1.5)
xlabel('p'); ylabel('SSdU'); grid on
legend('Data','Cheb5 (fminunc)','Location','southeast')

figure;
plot(f, SSE, 'LineWidth', 1.5); hold on
plot(f, SSE_fit, '--', 'LineWidth', 1.5)
xlabel('p'); ylabel('SSE'); grid on
legend('Data','Cheb5 (fminunc)','Location','southeast')

figure;
plot(f, time, 'LineWidth', 1.5); hold on
plot(f, time_fit, '--', 'LineWidth', 1.5)
xlabel('p'); ylabel('Time'); grid on
legend('Data','Cheb2 constrained (fminunc)','Location','best')

% -----------------------------
% Coefficients
% -----------------------------
fprintf('Cheb5 coeffs SSdU (c0..c5):\n');
fprintf('% .6e  ', c_SSdU);
fprintf('\n\n');

fprintf('Cheb5 coeffs SSE  (c0..c5):\n');
fprintf('% .6e  ', c_SSE);
fprintf('\n');


fprintf('Cheb2 constrained time coeffs:\n');
fprintf('% .6e  ', c_time);
fprintf('\n');

%% FUNCTIONS
Cheb5([-1 0 1], c_SSdU)
Cheb5([-1 0 1], c_SSE)
Cheb_time([-1 0 1], c_time)

% =============================
% Objectives
% =============================
function J = obj_Cheb5(c, x, y)
    r = y - Cheb5(x, c);
    J = r.'*r;
end

function J = obj_Cheb2_constrained(c, x, y)
    r = y - Cheb_time(x, c);
    J = r.'*r;
end

% =============================
% Models
% =============================
function y = Cheb5(x, c)
    c = c(:);
    if numel(c) ~= 6
        error('Cheb5 expects c = [c0;c1;c2;c3;c4;c5].')
    end

    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;
    T4 = 8*x.^4 - 8*x.^2 + 1;
    T5 = 16*x.^5 - 20*x.^3 + 5*x;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3 + c(5)*T4 + c(6)*T5;
end


function y = Cheb_time(x, c)
    c = c(:);

    c0 = c(1);
    c1 = c(2);
    c2 = c(3);
    c3 = c(4);

    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;

    y = c0*T0 + c1*T1 + c2*T2 + c3*T3;


    k = 20;                                  % smoothness (10–50 typical)
    softplus = @(z) log1p(exp(k*z)) / k;     % stable for moderate k
    y = y + softplus(-y);                    % identity when f>=0, ~0 when f<0

end