function  [dx] = dilution_reduced(t,x,u,par)
    % inputs:
    % t=time, 
    % x: state vector 
    % u: input vector (airflow, Fin, Fout, ..)
    % par: struct for the parameters of the model, can be fed as
    % Filter.augm.par if we need to change parameters in the current model's
    % call

    % outputs: the left side of the ODEs, i.e. the dxdt

    % States
    V = x(1); X = x(2); S = x(3);
    % Manipulated variables
    u = u(:)';
    Fin = u(1:end-1)';
    Fout = u(end);

    Sin = zeros(1,length(Fin));
    Sin(1)= par.Sin;

    % Parameters
    Y_XSinv   = par.Y_XSinv;
    mu_max = par.mu;   % h^-1
    Ks    = par.Ks;   % g L^-1 
    kd     = par.kd;

    mu = mu_max.*(S ./(Ks + S));   % h^-1

    % Differential equations:
    dV   = sum(Fin) - Fout;
    dX   = ( -(sum(Fin)/V) + mu  - kd ).*X ; %Biomass
    dS   = (Sin - S)*(Fin/V) - mu .*X .*Y_XSinv ; % Substrate

    % Initial empty vector
    dx = x*0; % This creates an empty output column vector
    
    % Output:
    dx(1)= dV; 
    dx(2)= dX; 
    dx(3)= dS;
end