syms Fin Fd Fout X V S

syms Y_XSinv mu_max Ks kd Sin

    mu = mu_max.*(S ./(Ks + S));   % h^-1

    % Differential equations:
    dV   = Fin + Fd - Fout;
    dX   = ( -(Fin+Fd)/V + mu  - kd ).*X ; %Biomass
    dS   = (Sin - S)*(Fin/V) -S*Fd/V - mu .*X .*Y_XSinv ; % Substrate

    % Initial empty vector    
    % Output:
    dXdt(1)= dV; 
    dXdt(2)= dX; 
    dXdt(3)= dS;


X = [V X S];
U = [Fin; Fd; Fout];

    
    % Symbolic state space
    Asym = jacobian(dXdt,X)
    Bsym = jacobian(dXdt,U)