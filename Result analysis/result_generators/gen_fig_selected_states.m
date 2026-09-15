function gen_fig_selected_states(ctx)
%GEN_FIG_SELECTED_STATES State trajectories of the selected controllers.
%   Result: results/graphical_results/selected_states.png/.pdf
%
%   Diagnostic grid, 3 states by one column per controller, columns ordered by
%   increasing J_track (see ra_selected_controllers): the best single-fidelity
%   point, the 3 lowest-J_track multi-fidelity points, and every damped-Rdu
%   blend in results/rdu_damping/. Scenario 1 solid,
%   scenario 2 dashed, setpoint grey dotted. Every trace runs to its own
%   i_last, so the multi-fidelity columns stop short of 10 h; the x limit is
%   pinned at the full horizon so that is visible rather than rescaled away.
%
%   Not a paper figure. It exists to be looked at before deciding what to
%   merge into one.

Sel        = ra_selected_controllers(ctx);
n          = numel(Sel);
stateNames = ["V", "X", "S"];
stateUnits = ["L", "g/L", "g/L"];
tfH        = 10;

fig = figure('Color', 'w', 'Name', 'Selected controllers - states');
tiledlayout(fig, 3, n, 'Padding', 'compact', 'TileSpacing', 'compact');

for r = 1:3
    for k = 1:n
        ax = nexttile((r - 1) * n + k);
        hold(ax, 'on');
        o = Sel(k).out;
        T = double(o.T(:));
        col = ra_selected_color(ctx, Sel(k).arm);

        for s = 1:numel(o.case)
            c = o.case(s);
            if s == 1; ls = '-'; else; ls = '--'; end
            plot(ax, T, double(c.Ysp(:, r)), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
            plot(ax, T, double(c.Y(:, r)), ls, 'Color', col, 'LineWidth', 1.4);
        end

        xlim(ax, [0 tfH]);
        if r == 1
            % the latex interpreter is the groot default, so the underscore in
            % "BO_1" has to be escaped or it is read as a subscript command.
            title(ax, sprintf('%s  ($z=%.3f$)', ra_tex_label(Sel(k).label), Sel(k).z), ...
                'Interpreter', 'latex');
        end
        if k == 1
            ylabel(ax, sprintf('%s (%s)', stateNames(r), stateUnits(r)));
        end
        if r == 3
            xlabel(ax, '$t$ (h)');
        end
        box(ax, 'off');
        set(ax, 'FontSize', 11);
    end
end

sgtitle(fig, 'Selected controllers, states. Solid: scenario 1. Dashed: scenario 2. Grey: setpoint.', ...
    'Interpreter', 'latex', 'FontSize', 12);

save_plot_outputs(fig, fullfile(ctx.graphicsDir, 'selected_states.png'), ...
    11, 320 * n, 720);
end
