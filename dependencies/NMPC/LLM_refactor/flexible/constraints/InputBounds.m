classdef InputBounds < handle
    % InputBounds  Hard bounds on the `u` block.

    properties (SetAccess = private)
        umin
        umax
    end

    methods
        function obj = InputBounds(umin, umax)
            arguments
                umin (1,:) double
                umax (1,:) double
            end
            assert(numel(umin) == numel(umax), ...
                'InputBounds:size', 'umin/umax size mismatch.');
            obj.umin = umin;
            obj.umax = umax;
        end

        function ov = bound_overrides(obj, layout)
            sz = layout.size_of('u');   % [m, nu]
            ov.u = struct( ...
                'lo', repmat(obj.umin, sz(1), 1), ...
                'hi', repmat(obj.umax, sz(1), 1));
        end
    end
end
