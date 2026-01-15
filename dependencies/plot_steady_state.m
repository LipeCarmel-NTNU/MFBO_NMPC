function plot_steady_state(par)
    % Generation
    linewidth = 1.5;
    N = 2000; % Resolution
    X_start = 0;
    X_end = 30;

    X = linspace(X_start, X_end, N);

    % Finer resolution at the intercept
    fine_start = 28.5; fine_end = 28.9;
    X = [linspace(0, fine_start, floor(N/4)),...
        linspace(fine_start, fine_end, floor(N/2)),...
        linspace(fine_end, X_end, floor(N/4))];

    N = length(X);

    vS = zeros(N, 2);
    for i = 1 : N
        [S_low, S_high] = steady_state(par, X(i));
        vS(i, 1) = S_low;
        vS(i, 2) = S_high;
    end

    %% Plotting

    figure
    plot(X, vS(:, 1), 'b--', 'LineWidth',linewidth)
    hold on
    plot(X, vS(:, 2), 'r-', 'LineWidth',linewidth)
    set(gca, 'YScale', 'log')
    
    xlabel('Biomass (g/L)', 'Interpreter','latex')
    ylabel('Susbtrate (g/L)', 'Interpreter','latex')
    set_font_size()

    % lgd = legend('Lower subtrate steady-state', 'Higher subtrate steady-state');
    % lgd.FontSize = 16;
    ylim([0 inf])
    xlim([X_start X_end])

    set_fig_size()


    function [S_low, S_high] = steady_state(par, X)
        %% Parameters
        % Constants
        V = 1; % so Fin = D
        Sin = par.Sin;

        % Model parameters
        Y_XSinv   = par.Y_XSinv;
        mu_max = par.mu_max;
        Ks    = par.Ks;   % g L^-1
        Y_CO2Xinv = 1/par.Y_CO2X;
        kd     = par.kd;


        %% Nullclines
        S1 = solve_starvation(X);
        S2 = solve_abundance(X);

        % Pick the correct solution
        % % (they may switch in near-maximum X conditions)
        S_low = min([S1 S2]);
        S_high = max([S1 S2]);

        function S_sol = solve_starvation(X)
            S = kd * Ks / (mu_max - kd); % D = 0 -> mu = kd
            D = kd * X * Y_XSinv / (Sin - S);

            x0 = [S, D];
            [x, ~, flag] = fsolve(@(x) ode_system(X, x(1), x(2)), x0, optimset('Display','off', 'MaxFunEvals', 1e4));
            if flag < 1
                % Not solved
                S_sol = [];
            else
                S_sol = x(1);
            end
        end
        function S_sol = solve_abundance(X)
            D = mu_max - kd;
            S = Sin - (mu_max * X * Y_XSinv)/ D;

            x0 = [S, D];
            [x, ~, flag] = fsolve(@(x) ode_system(X, x(1), x(2)), x0, optimset('Display','off', 'MaxFunEvals', 1e4));
            if flag < 1
                % Not solved
                S_sol = [];
            else
                S_sol = x(1);
            end
        end
        function dxdt = ode_system(X, S, Fin)
            dX   = -X.*(Fin/V) + mu_max .*(S ./(Ks + S)) .*X - kd .*X ;
            dS   = (Sin-S).*(Fin/V) - mu_max *(S ./(Ks + S)) .*X .*Y_XSinv;
            dxdt = [dX, dS];
        end
    end

end