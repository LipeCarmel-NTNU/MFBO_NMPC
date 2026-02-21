function colors = nature_methods_colors(n)
%NATURE_METHODS_COLORS Wong (2011) Nat Methods colorblind-safe palette.
%   If called with no inputs, returns a struct of named colors (RGB in [0,1]).
%   Otherwise returns an n-by-3 matrix in this sequence:
%     blue, bluish green, reddish purple, vermillion, orange, sky blue, yellow
%
% Source: Wong, B. Points of view: Color blindness. Nat Methods 8, 441 (2011).
% https://doi.org/10.1038/nmeth.1618
% 
% See Fig. 2: https://www.nature.com/articles/nmeth.1618/figures/2

    % ---- Define palette as a struct (paper order; RGB in [0,1]) ----
    % NATURE_COLOR.Black         = [  0,   0,   0] / 255;  % omitted
    NATURE_COLOR.Orange        = [230, 159,   0] / 255;
    NATURE_COLOR.SkyBlue       = [ 86, 180, 233] / 255;
    NATURE_COLOR.BluishGreen   = [  0, 158, 115] / 255;
    NATURE_COLOR.Yellow        = [240, 228,  66] / 255;
    NATURE_COLOR.Blue          = [  0, 114, 178] / 255;
    NATURE_COLOR.Vermillion    = [213,  94,   0] / 255;
    NATURE_COLOR.ReddishPurple = [204, 121, 167] / 255;

    % If no inputs, return the struct
    if nargin == 0
        colors = NATURE_COLOR;
        return;
    end

    % Validate: must be an integer scalar and not exceed the order length
    if ~isscalar(n) || n ~= floor(n)
        error('n must be an integer scalar.');
    end

    order = {'Blue', 'BluishGreen', 'ReddishPurple', 'Vermillion', ...
             'Orange', 'SkyBlue', 'Yellow'};

    if n > numel(order)
        error('n must be <= %d.', numel(order));
    end

    % Build output matrix in the requested sequence
    colors = zeros(n, 3);
    for i = 1:n
        colors(i, :) = NATURE_COLOR.(order{i});
    end
end
