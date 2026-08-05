% Alternative to the surrogate model: start with an arbitrary polynomial
% and transform to generate the necessary properties.
% Issue: biased towards y=x because the slope is derived from integration.
% - flat: final curve is 2x - x^2
% - unit step: slightly more than y=x, flattening in the end
% - high then low: quick saturation
% i.e: the integral needs to be concentrated near the start to result in
% appreciable difference
close all

x = linspace(0, 1, 100);


% Evaluate the polynomial at the points in x: coeff applied to shifted
% Legendre polynomials P_n(2x - 1), orthogonal on [0, 1] (uniform weight),
% instead of the trivial monomial basis.
coeff = 10*randn([1, 3]);
% t = 2*x - 1;                      % map [0, 1] -> [-1, 1]
% P = zeros(numel(coeff), numel(x));
% P(1, :) = 1;                      % P_0
% if numel(coeff) > 1
%     P(2, :) = t;                  % P_1
% end
% for n = 2:numel(coeff) - 1        % Bonnet recurrence: row n+1 holds P_n
%     P(n+1, :) = ((2*n - 1) * t .* P(n, :) - (n - 1) * P(n-1, :)) / n;
% end
% y = coeff * P;

y = x;
% y = 2*ones(size(x));
% y(round(length(y)/4): end) = 0; % heavy on the start
%y(1:round(length(y)/2)) = 0;      % heavy on the end

figure
plot(x,y)
hold on
yline(0)

sqry = y.^2;
figure
plot(x,sqry)
hold on
yline(0)

I = cumtrapz(x, sqry);
I = 1 - I/I(end);
figure
plot(x,I)
hold on
yline(0)


f = cumtrapz(I);
f = f/f(end);
figure
plot(x,f)
hold on
yline(0)
plot([0 1], [0 1], 'k-')
