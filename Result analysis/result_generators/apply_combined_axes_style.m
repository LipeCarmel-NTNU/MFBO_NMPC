function apply_combined_axes_style(ax, fontSize, xLim, yLim)
%APPLY_COMBINED_AXES_STYLE Shared log-log axes/tick styling for combined views.
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
xlim(ax, xLim);
ylim(ax, yLim);
xlabel(ax, '$J_{\mathrm{TV}}$');
ylabel(ax, '$J_{\mathrm{track}}$');
grid(ax, 'off');
box(ax, 'off');
axes(ax);
format_tick(1, 1);
end
