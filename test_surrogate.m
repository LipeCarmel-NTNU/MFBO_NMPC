% compare_models_to_data.m
close all; clc

load('surrogate_data_2.mat')

% -----------------------------
% Rebuild p and data
% -----------------------------
N = length(out.case(1).SSE);
f = (1:N).' / N;
x = 2*f - 1;

m = out.theta(3) + 1;
p = out.theta(2) +  m;

time = cumsum(out.case(1).RUNTIME) + cumsum(out.case(2).RUNTIME);

SSdU = [0; cumsum(out.case(1).SSdU)] + [0; cumsum(out.case(2).SSdU)];
SSdU = SSdU / SSdU(end);

SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
SSE  = SSE / SSE(end);

% -----------------------------
% Model coefficients
% -----------------------------
c_SSdU = [ ...
    6.348568e-01
    4.707209e-01
   -1.352747e-01
    2.693132e-02
    1.620533e-02
   -9.064012e-03 ];

c_SSE = [ ...
    7.833013e-01
    3.766002e-01
   -2.445100e-01
    1.103409e-01
   -2.894632e-02
    5.486239e-04 ];

c_time = [3.641170e+03   4.295613e+03   6.326895e+02   2.022585e+01];

% -----------------------------
% Model evaluation
% -----------------------------
SSdU_f = Cheb5(x, c_SSdU);
SSE_f  = Cheb5(x, c_SSE);
time_hat = max(Cheb2(x, c_time), 1);

SSdU_f = min(max(SSdU_f, 0), 1);
SSE_f  = min(max(SSE_f,  0), 1);

% -----------------------------
% Plots: data vs model
% -----------------------------
figure;
plot(f, SSdU, 'LineWidth', 1.5); hold on
plot(f, SSdU_f, '--', 'LineWidth', 1.5)
xlabel('f'); ylabel('SSdU'); grid on
legend('Data','Model','Location','southeast')

figure;
plot(f, SSE, 'LineWidth', 1.5); hold on
plot(f, SSE_f, '--', 'LineWidth', 1.5)
xlabel('f'); ylabel('SSE'); grid on
legend('Data','Model','Location','southeast')

figure;
plot(f, time, 'LineWidth', 1.5); hold on
plot(f, time_hat, '--', 'LineWidth', 1.5)
xlabel('f'); ylabel('Time'); grid on
legend('Data','Model','Location','best')

figure; plot(f, time./time_hat)

% =============================
% Models
% =============================
function y = Cheb5(x, c)
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;
    T4 = 8*x.^4 - 8*x.^2 + 1;
    T5 = 16*x.^5 - 20*x.^3 + 5*x;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3 + c(5)*T4 + c(6)*T5;
end

function y = Cheb2(x, c)
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2;
end