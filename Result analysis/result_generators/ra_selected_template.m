function S = ra_selected_template()
%RA_SELECTED_TEMPLATE One empty entry of the selected-controller struct array.
%   ra_selected_controllers and ra_rdu_controllers both build their entries
%   from this, so the two arrays concatenate: MATLAB requires identical field
%   names in identical order.
S = struct('label', "", 'arm', "", 'caseName', "", 'file', "", ...
           'out', [], 'z', NaN, 'SSE', NaN, 'SSdU', NaN, ...
           'Np', NaN, 'Nc', NaN);
end
