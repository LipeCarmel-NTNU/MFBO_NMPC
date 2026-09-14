function gen_fig_selected_inputs(ctx)
%GEN_FIG_SELECTED_INPUTS Input trajectories of the selected controllers.
%   Result: results/graphical_results/selected_inputs.png/.pdf
%
%   Same grid as gen_fig_selected_states: 3 inputs by one column per
%   controller, columns ordered by increasing J_track. u1 and u2 are the feed
%   flows, u3 the outflow. Drawn with stairs because the input is held over
%   the sample. Scenario 1 solid, scenario 2 dashed.
%
%   This is where the J_TV gap lives: the multi-fidelity columns move the
%   inputs an order of magnitude more than the baseline ones do.

Sel        = ra_selected_controllers(ctx);
n          = numel(Sel);
inputNames = ["u_1 = F_{in,1}", "u_2 = F_{in,2}", "u_3 = F_{out}"];
tfH        = 10;

uMax = 0;
for k = 1:n
    for s = 1:numel(Sel(k).out.case)
        uMax = max(uMax, max(double(Sel(k).out.case(s).U), [], 'all'));
    end
end
uMax = max(uMax, eps);

fig = figure('Color', 'w', 'Name', 'Selected controllers - inputs');
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
            stairs(ax, T, double(c.U(:, r)), ls, 'Color', col, 'LineWidth', 1.2);
        end

        xlim(ax, [0 tfH]);
        ylim(ax, [0 1.05 * uMax]);
        if r == 1
            title(ax, sprintf('%s  ($N_c=%d$)', ra_tex_label(Sel(k).label), Sel(k).Nc), ...
                'Interpreter', 'latex');
        end
        if k == 1
            ylabel(ax, "$" + inputNames(r) + "$ (L/h)");
        end
        if r == 3
            xlabel(ax, '$t$ (h)');
        end
        box(ax, 'off');
        set(ax, 'FontSize', 11);
    end
end

sgtitle(fig, 'Selected controllers, inputs. Solid: scenario 1. Dashed: scenario 2.', ...
    'Interpreter', 'latex', 'FontSize', 12);

save_plot_outputs(fig, fullfile(ctx.graphicsDir, 'selected_inputs.png'), ...
    11, 320 * n, 720);
end
