function varargout = linearize(xss, uss, model)
    % Linearize time invariant ODE
    arguments
        xss double
        uss double
        model function_handle
    end

    nx = length(xss);
    nu = length(uss);

    % Symbolic state space
    [Asym, Bsym, X, U] = symbolic_ss(model, nx, nu);
    
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