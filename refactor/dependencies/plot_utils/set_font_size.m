function set_font_size(newFontSize)
    arguments
        newFontSize = 14
    end

    grid off
    box off
    % Get the current figure
    fig = gcf;

    % Find all axes objects in the figure
    ax = findall(fig, 'type', 'axes');

    % Update the font size for all axes and text elements
    for i = 1:length(ax)
        ax(i).FontSize = newFontSize;
        % If you want to specifically set axes titles, labels, etc.:
        ax(i).Title.FontSize = newFontSize;
        ax(i).XLabel.FontSize = newFontSize;
        ax(i).YLabel.FontSize = newFontSize;
        ax(i).ZLabel.FontSize = newFontSize; % If there's a z-axis
    end

end