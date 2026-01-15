function colors = good_colors(n)
%GOOD_COLORS Return an n-by-3 matrix of distinct RGB colors.
%   The base order is: red, green, blue, purple, yellow, orange.
%   Colors are dropped in this sequence until n remain:
%   green -> orange -> purple -> yellow.
%
%   n must be an integer in [1, 6].

    arguments
        n int64
    end

    if ~isscalar(n) || n ~= floor(n)
        error('n must be an integer scalar.');
    end
    if n < 1 || n > 6
        error('n must be between 2 and 6 to respect the specified drop sequence.');
    end

    % Define colors in RGB (0–255)
    % https://www.molecularecologist.com/2020/04/23/simple-tools-for-mastering-color-in-scientific-figures/
    C.red    = [255,  31,  91];  % #FF1F5B
    C.green  = [  0, 205, 108];  % #00CD6C
    C.blue   = [  0, 154, 222];  % #009ADE
    C.purple = [175,  88, 186];  % #AF58BA
    C.yellow = [255, 198,  30];  % #FFC61E
    C.orange = [242, 133,  34];  % #F28522
    

    % Base order to keep in the output
    order = {'red','green','blue','purple','yellow','orange'};

    % Drop sequence as specified
    drop_seq = {'green','orange','purple','yellow', 'red'};

    % Number of colors to drop
    kdrop = numel(order) - n;

    % Remove colors according to drop sequence, preserving order
    for i = 1:kdrop
        name = drop_seq{i};
        idx = find(strcmp(order, name), 1, 'first');
        if ~isempty(idx)
            order(idx) = [];
        end
    end

    % Build output matrix in remaining order and scale to [0,1]
    colors = zeros(n, 3);
    for i = 1:n
        colors(i, :) = C.(order{i}) / 255;
    end
end
