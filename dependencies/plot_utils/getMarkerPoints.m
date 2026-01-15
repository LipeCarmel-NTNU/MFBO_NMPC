function [marker_t, marker_y] = getMarkerPoints(t, y, n_markers, jitter_percent)
    % getMarkerPoints - Select marker points with random jitter
    %
    % Inputs:
    %   t              - Time or x-axis vector
    %   y              - Data or y-axis vector
    %   n_markers      - Number of markers to place (default: 5)
    %   jitter_percent - Random horizontal jitter as % of spacing (default: 50)
    % Outputs:
    %   marker_t - Time/x values for markers
    %   marker_y - Data/y values for markers
    %
    % Example:
    %   % Create two sinusoidal functions
    %   t = linspace(0, 4*pi, 200);
    %   y1 = sin(t);
    %   y2 = cos(t);
    %
    %   % Plot with markers
    %   figure;
    %   plot(t, y1, 'b-', 'LineWidth', 1.5); hold on;
    %   plot(t, y2, 'r-', 'LineWidth', 1.5);
    %
    %   % Add markers with jitter
    %   [mt1, my1] = getMarkerPoints(t, y1, 6, 40);
    %   [mt2, my2] = getMarkerPoints(t, y2, 6, 40);
    %
    %   plot(mt1, my1, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    %   plot(mt2, my2, 'rs', 'MarkerSize', 8, 'LineWidth', 2);
    %
    %   legend('sin(t)', 'cos(t)', 'sin markers', 'cos markers');
    %   grid on;
    
    arguments
        t double
        y double
        n_markers (1,1) {mustBeInteger, mustBePositive} = 5
        jitter_percent (1,1) double {mustBeNonnegative} = 50
    end
    
    % Get time range
    t_min = min(t);
    t_max = max(t);
    
    % Generate evenly spaced time points
    % linspace gives exactly n_markers points from t_min to t_max
    t_calc = linspace(t_min, t_max, n_markers);
    
    % Calculate time step for jitter calculation
    if n_markers > 1
        t_step = (t_max - t_min) / (n_markers - 1);
    else
        t_step = (t_max - t_min) / 2;  % For single marker, use half range
    end
    
    % Add random jitter proportional to time spacing
    jitter = (rand(size(t_calc)) - 0.5) * 2 * t_step * jitter_percent / 100;
    t_calc = t_calc + jitter;
    
    % Ensure time points remain within valid range
    t_calc = max(t_min, min(t_max, t_calc));
    
    % Ensure uniqueness and sorted order
    t_calc = unique(sort(t_calc));
    
    % Find closest indices to the calculated time points
    idx = zeros(size(t_calc));
    for i = 1:length(t_calc)
        [~, idx(i)] = min(abs(t - t_calc(i)));
    end
    
    % Ensure uniqueness of indices
    idx = unique(idx);
    
    % Extract marker points
    marker_t = t(idx);
    marker_y = y(idx);
end