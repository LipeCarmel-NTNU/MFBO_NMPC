function [Asym, Bsym, X, U] = symbolic_ss(model, nx, nu)
    arguments
        model function_handle
        nx double {mustBeNonnegative,mustBeInteger}
        nu double {mustBeNonnegative,mustBeInteger}
    end

    X = sym('X',[nx 1]);
    U = sym('U',[nu 1]);

    % Symbolic ODE
    if nu > 0
        dXdt = model(X, U);
    else
        dXdt = model(X);
    end

    Asym = jacobian(dXdt, X);
    Bsym = jacobian(dXdt, U);

end