function save_figs(fig_group_number, name, varargin)
    % Parse optional file extension
    if nargin > 2 && ischar(varargin{1})
        ext = varargin{1};
        if ~startsWith(ext, '.')
            ext = ['.', ext]; % Ensure leading dot
        end
    else
        ext = '.pdf'; % Default extension
    end

    fig_number = fig_group_number;
    i = 1;

    while ishandle(fig_number)

        h = figure(fig_number);
        h.InnerPosition = [2900 400 560 420];
        h.OuterPosition = [2900 400 580 515];

        filename = strcat(name, string(i), ext);
        save_figure(fig_number, filename)

        i = i + 1;
        fig_number = fig_number + 1;
    end
end