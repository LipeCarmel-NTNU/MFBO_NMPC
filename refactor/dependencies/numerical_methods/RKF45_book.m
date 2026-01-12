function [t,y] = RKF45_book(system, tspan, x0, h)
% NOTE: THE METHOD WAS MODIFIED TO PREVENT NEGATIVE VALUES. THIS FUNCTION
% CANNOT BE USED FOR SOLUTIONS WITH NEGATIVE VALUES.
%
% Runge-Kuta-Fehlberg Method adapted from page 466 of:
%
% Numerical Methods Using MATLAB (3rd Edition)
% Mathews, John H., Fink, Kurtis D.
% Published by Prentice Hall College Div (1998)
% ISBN 10: 0132700425  ISBN 13: 9780132700429
    if size(x0, 2) == 1
        x0 = x0';
    end
    t0 = tspan(1);
    tf = tspan(2);
    
    % Coefficients for RK45
    a = [0 1/4 3/8 12/13 1 1/2];
    b = [0 0 0 0 0;
         1/4 0 0 0 0;
         3/32 9/32 0 0 0;
         1932/2197 -7200/2197 7296/2197 0 0;
         439/216 -8 3680/513 -845/4104 0;
         -8/27 2 -3544/2565 1859/4104 -11/40];
    c = [16/135 0 6656/12825 28561/56430 -9/50 2/55];
    c_star = [25/216 0 1408/2565 2197/4104 -1/5 0];
    
    t = t0;
    y = x0';
    
    tol = 1e-8;
    safety = 0.9;
    
    T = t0;
    Y = x0;
    
    if nargin > 3
        h_new = h;
    else
        h_new = 0.1 * norm(x0) / norm(system(t0, x0'));
    end
    
    while t < tf
        h = h_new;
        if t + h > tf
            h = tf - t;
        end

        k1 = h * system(t, y);
        k2 = h * system(t + a(2)*h, y + b(2,1)*k1);
        k3 = h * system(t + a(3)*h, y + b(3,1)*k1 + b(3,2)*k2);
        k4 = h * system(t + a(4)*h, y + b(4,1)*k1 + b(4,2)*k2 + b(4,3)*k3);
        k5 = h * system(t + a(5)*h, y + b(5,1)*k1 + b(5,2)*k2 + b(5,3)*k3 + b(5,4)*k4);
        k6 = h * system(t + a(6)*h, y + b(6,1)*k1 + b(6,2)*k2 + b(6,3)*k3 + b(6,4)*k4 + b(6,5)*k5);
        
        y_new =  y + c(1)*k1 + c(2)*k2 + c(3)*k3 + c(4)*k4 + c(5)*k5 + c(6)*k6;
        y_star = y + c_star(1)*k1 + c_star(2)*k2 + c_star(3)*k3 + c_star(4)*k4 + c_star(5)*k5 + c_star(6)*k6;
        
        % Error estimate (maximum absolute error)
        %e = max(norm(y_new - y_star, inf), - min([y_new(y_new<0); 0]));
        e = norm(y_new - y_star, inf);

        % Step size control
        h_new = h * safety * 0.84*(tol / e)^(1/4);

    % Keep integrating if whitin tol
        if e < tol
            y_new(y_new<0) = 0; % ENSURING POSITIVE VALUES      
            t = t + h;
            y = y_new;
            T = [T; t];
            Y = [Y; y'];
        end
    end
    t = T;
    y = Y;
end