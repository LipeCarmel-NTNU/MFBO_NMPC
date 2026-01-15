function set_fig_size(W, L)
    arguments
        W = 650
        L = 400
    end
    % Used to generate figures with the same size
    fig = gcf;
    set(fig, 'Units', 'pixels');
    fig_size = get(fig, 'Position');
    set(fig, 'Position', [fig_size(1) fig_size(2) W L]);
end