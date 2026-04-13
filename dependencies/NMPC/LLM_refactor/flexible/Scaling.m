classdef Scaling < handle
    % Scaling  Per-block multiplicative scale factors.
    %
    %   Stores a positive row vector per registered block. Unregistered
    %   blocks default to ones (identity). Conversion convention:
    %       w_scaled = x_phys ./ scale        (pack)
    %       x_phys   = w_scaled .* scale      (unpack)

    properties (Access = private)
        factors = struct()    % block name -> 1×n_cols positive row vector
    end

    methods
        function set(obj, name, factor)
            name = char(name);
            factor = reshape(factor, 1, []);
            assert(all(factor > 0) && all(isfinite(factor)), ...
                'Scaling:pos', 'Scale factors must be finite and positive.');
            obj.factors.(name) = factor;
        end

        function f = get(obj, name, n_cols)
            % Return a 1×n_cols row. If not registered, returns ones(1, n_cols).
            name = char(name);
            if isfield(obj.factors, name)
                f = obj.factors.(name);
                if nargin >= 3
                    assert(numel(f) == n_cols, 'Scaling:size', ...
                        'Scale for "%s" has %d entries, expected %d.', ...
                        name, numel(f), n_cols);
                end
            else
                if nargin < 3
                    n_cols = 1;
                end
                f = ones(1, n_cols);
            end
        end

        function tf = has(obj, name)
            tf = isfield(obj.factors, char(name));
        end
    end
end
