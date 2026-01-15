function  [dx] = dilution_Monod_model(t,x,u,par,adaptPar)
    % inputs:
    % t=time, 
    % x: state vector 
    % u: input vector (airflow, Fin, Fout, ..)
    % par: struct for the parameters of the model, can be fed as
    % Filter.augm.par if we need to change parameters in the current model's
    % call
 if nargin>4 && adaptPar
     % the adaptive parameters
    alpha=par.ModAdapt.mu;
    beta=par.ModAdapt.beta;
 else
     alpha=1;beta=1; % no affect
 end

    % outputs: the left side of the ODEs, i.e. the dxdt

    % States
    V = x(1); X = x(2); S = x(3); CO2 = x(4);
    % Manipulated variables
    u = u(:)';
    Q = u(1);
    Fin = u(2:end-1)';
    Fout = u(end);

    Sin = zeros(1,length(Fin));
    Sin(1)= par.Sin;

    % Parameters
    Y_XSinv   = par.Y_XSinv;
    mu = par.mu;   % h^-1
    mu_max = par.mu_max;
    Ks    = par.Ks;   % g L^-1 
    Y_CO2Xinv = par.Y_CO2Xinv;
    kd     = par.kd;

    % Differential equations:
    dV   = sum(Fin) - Fout;
    dX   = -X.*(sum(Fin)/V) + alpha*mu .*(S ./(Ks + S)) .*X - kd .*X ; %Biomass
    dS   = (Sin-S)*(Fin/V) - mu *(S ./(Ks + S)).*X .*Y_XSinv ; % Substrate
    dCO2 = -CO2.*Q + (1 .*(S ./(Ks + S)) .*X).*Y_CO2Xinv *beta; %CO2

    % Initial empty vector
    dx = x*0; % This creates an empty output column vector
    
    % Output:
    dx(1)= dV; 
    dx(2)= dX; 
    dx(3)= dS;
    dx(4)= dCO2;
end