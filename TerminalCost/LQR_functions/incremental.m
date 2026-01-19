function [Ai, Bi] = incremental(A,B,Ts)
    % Assumes C = I, D = 0;
    
    nx = size(A,1);
    nu = size(B,2);

    sysc = ss(A,B,eye(nx),zeros(nx,nu));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    Ai = [Ad, Bd;
          zeros(nu,nx), eye(nu)];
    Bi = [Bd;
          eye(nu)];
end