
clc
N = 100;

% Case 1
out.case(1).RUNTIME = 0.5 + 0.1*rand(N,1);
out.case(1).SSdU    = abs(randn(N-1,1));
out.case(1).SSE     = abs(randn(N,1));

% Case 2
out.case(2).RUNTIME = 0.4 + 0.1*rand(N,1);
out.case(2).SSdU    = abs(randn(N-1,1));
out.case(2).SSE     = abs(randn(N,1));

% Cumulative quantities, summed over case 1 and 2
N = length(out.case(1).SSE);
p = (1:N).' / N;

time = cumsum(out.case(1).RUNTIME) + cumsum(out.case(2).RUNTIME)./p;

SSdU = [0; cumsum(out.case(1).SSdU)] + ...
       [0; cumsum(out.case(2).SSdU)];
SSdU = SSdU./p;

SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
SSE = SSE./p;

% Normalised progress variable

% Plots
figure;
plot(p, SSdU, 'LineWidth', 1.5);
xlabel('p');
ylabel('SSdU');
grid on;

figure;
plot(p, SSE, 'LineWidth', 1.5);
xlabel('p');
ylabel('SSE');
grid on;

figure;
plot(p, time, 'LineWidth', 1.5);
xlabel('p');
ylabel('Time');
grid on;
