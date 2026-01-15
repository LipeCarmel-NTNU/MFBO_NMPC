function varargout = linearize(xss, uss, model)
    % Linearize time invariant ODE
    arguments
        xss double
        uss double
        model function_handle
    end

    nx = length(xss);
    nu = length(uss);
    
    X = sym('X',[1 nx]).';
    U = sym('U',[1 nu]).';
    
    % Symbolic ODE
    if nu > 0
        dXdt = model(X, U);
    else
        dXdt = model(X);
    end
    
    % Symbolic state space
    Asym = jacobian(dXdt,X);
    Bsym = jacobian(dXdt,U);
    
    % State space
    A = double(subs(Asym, [X; U], [xss'; uss']));
    B = double(subs(Bsym, [X; U], [xss'; uss']));

    if nargout == 2
        varargout = {A, B};
        return 
    end

    % Laplace transform
    s = tf('s');
    G = (s * eye(size(A,1)) - A) \ B;

    if nargout == 1
        varargout = G;
    else
        varargout = {A, B, G};
    end
end