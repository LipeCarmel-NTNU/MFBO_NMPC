function plot_simulation(out, case_id)
    %PLOT_SIMULATION Plot states and inputs for a given simulation case.
    % Utility function (not called by the main script). Use after loading a saved out_*.mat.

    Y   = out.case(case_id).Y;
    Ysp = out.case(case_id).Ysp;
    U   = out.case(case_id).U;
    dt  = out.case(case_id).dt;
 
    N = size(Y,1);
    T = 0 : dt : (N - 1)*dt;
    TF = 10;

    figure(10*case_id + 1);
    clf

    subplot(3,1,1);
    plot(T(1:N), Y(1:N,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)', 'Interpreter','latex');
    ylabel('State 1', 'Interpreter','latex');
    legend('Location','best', 'Interpreter','latex');
    xlim([0 TF])
    hold off;

    subplot(3,1,2);
    plot(T(1:N), Y(1:N,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)', 'Interpreter','latex');
    ylabel('State 2', 'Interpreter','latex');
    legend('Location','best', 'Interpreter','latex');
    xlim([0 TF])
    hold off;

    subplot(3,1,3);
    plot(T(1:N), Y(1:N,3), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,3), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)', 'Interpreter','latex');
    ylabel('State 3', 'Interpreter','latex');
    legend('Location','best', 'Interpreter','latex');
    xlim([0 TF])
    hold off;

    ax = findall(gcf, 'type', 'axes');
    for j = 1:numel(ax)
        ax(j).FontSize = 15;
        ax(j).XLabel.FontSize = 15;
        ax(j).YLabel.FontSize = 15;
    end

    figure(10*case_id + 2);
    clf
    plot(T(1:N), U(1:N,:), 'LineWidth', 2);
    grid on; box on;
    xlabel('Time (h)', 'Interpreter','latex');
    ylabel('Inputs', 'Interpreter','latex');
    xlim([0 TF])
end
