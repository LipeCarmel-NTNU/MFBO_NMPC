function save_figure(filename, fig_number, fig_bckp)
    %SAVE_FIGURE(fig_number, filename) Save figure without cutoff by adjusting paper size and export options
    %
    arguments
        filename = 'Figure.pdf'
        fig_number = NaN
        fig_bckp = true;
    end

    
    if isa(fig_number, 'matlab.ui.Figure')
        % Figure handle object passed directly
        fig = fig_number;
    elseif isnumeric(fig_number) && isscalar(fig_number) && ~isnan(fig_number) ...
            && ishandle(fig_number)
        % Numeric figure number
        fig = figure(fig_number);
    else
        % No valid target, fall back to current figure
        if ~(isnumeric(fig_number) && isscalar(fig_number) && isnan(fig_number))
            warning('Figure target not found; using current figure.');
        end
        fig = gcf;
    end

    % Ensure the figure uses its on-screen size when exporting
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperUnits', 'inches');

    % Adjust the paper size to match the figure
    fig_pos = get(fig, 'Position');            % [left bottom width height] in pixels
    dpi = get(0, 'ScreenPixelsPerInch');      % screen DPI
    widthIn = fig_pos(3) / dpi;
    heightIn = fig_pos(4) / dpi;
    set(fig, 'PaperSize', [widthIn, heightIn]);
    set(fig, 'PaperPosition', [0, 0, widthIn, heightIn]);
    try
        % Export using vector graphics without cutoff
        exportgraphics(fig, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    catch ME
        warning('FAILED TO SAVE %s: %s', filename, ME.message);
        % Fallback: use print with best-fit option
        print(fig, filename, '-bestfit', '-vector');
    end
    if fig_bckp
        [folder_path,name,extension]=fileparts(filename);
        fig_filename = fullfile(folder_path, [name, '.fig']);
        savefig(fig, fig_filename);
    end
end
