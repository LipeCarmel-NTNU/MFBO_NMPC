function format_tick(X_decimals, Y_decimals)

    ax = gca;
    hold on
    if ~isempty(X_decimals)
        ticks = ax.XTick;
        ax.XTickLabel = arrayfun(@(y) sprintf(['%.',num2str(X_decimals),'f'], y), ticks, 'UniformOutput', false);
    end
    
    if ~isempty(Y_decimals)
        ticks = ax.YTick;
        ax.YTickLabel = arrayfun(@(y) sprintf(['%.',num2str(Y_decimals),'f'], y), ticks, 'UniformOutput', false);
    end

    hold off
end