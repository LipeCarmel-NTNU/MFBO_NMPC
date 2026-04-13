classdef DecisionLayout < handle
    % DecisionLayout  Named-block layout of the decision vector w.
    %
    %   Each block is identified by a name and sized (n_steps x n_cols).
    %   Blocks are stored column-major inside w, in insertion order. After
    %   freeze(), pack/unpack are the only conversion API between physical
    %   named arrays (parts struct) and the scaled decision vector ws.
    %
    %   The layout itself is unitless; scaling is applied via a Scaling
    %   object supplied to pack/unpack.

    properties (SetAccess = private)
        names   = {}           % cell array of block names
        sizes   = zeros(0, 2)  % [n_steps, n_cols] per block
        offsets = []           % starting linear index (1-based) of each block
        total   = 0            % total length of w
        frozen  = false
    end

    methods
        function add(obj, name, n_steps, n_cols)
            % Register a new block. Only allowed before freeze().
            assert(~obj.frozen, 'DecisionLayout:frozen', ...
                'Layout is frozen; cannot add more blocks.');
            assert(ischar(name) || isstring(name), ...
                'DecisionLayout:name', 'Block name must be a string.');
            name = char(name);
            assert(~ismember(name, obj.names), ...
                'DecisionLayout:dup', 'Block "%s" already registered.', name);
            assert(n_steps >= 0 && n_cols >= 0, ...
                'DecisionLayout:size', 'Block sizes must be non-negative.');

            obj.names{end + 1}   = name;
            obj.sizes(end + 1,:) = [n_steps, n_cols];
            obj.offsets(end + 1) = obj.total + 1;
            obj.total            = obj.total + n_steps * n_cols;
        end

        function freeze(obj)
            obj.frozen = true;
        end

        function tf = has(obj, name)
            tf = ismember(char(name), obj.names);
        end

        function sz = size_of(obj, name)
            sz = obj.sizes(obj.index_of(name), :);
        end

        function [lo, hi] = range_of(obj, name)
            % Absolute 1-based index range of the block inside w.
            k  = obj.index_of(name);
            lo = obj.offsets(k);
            hi = lo + prod(obj.sizes(k,:)) - 1;
        end

        function names = list(obj)
            names = obj.names;
        end

        function n = total_size(obj)
            n = obj.total;
        end

        function w = pack(obj, parts, scaling)
            % parts : struct with one field per registered block
            %         (physical units, size n_steps×n_cols)
            % scaling : Scaling object (optional; identity if [])
            assert(obj.frozen, 'DecisionLayout:notfrozen', ...
                'Call freeze() before pack().');
            w = zeros(obj.total, 1);
            for k = 1 : numel(obj.names)
                nm   = obj.names{k};
                blk  = parts.(nm);
                sz   = obj.sizes(k, :);
                assert(isequal(size(blk), sz), ...
                    'DecisionLayout:packsize', ...
                    'Block "%s" expected size [%d %d], got [%d %d].', ...
                    nm, sz(1), sz(2), size(blk,1), size(blk,2));
                if ~isempty(scaling)
                    blk = blk ./ scaling.get(nm, sz(2));
                end
                w(obj.offsets(k) : obj.offsets(k) + prod(sz) - 1) = blk(:);
            end
        end

        function parts = unpack(obj, w, scaling)
            % Inverse of pack. Returns a struct of physical-unit arrays.
            assert(obj.frozen, 'DecisionLayout:notfrozen', ...
                'Call freeze() before unpack().');
            parts = struct();
            for k = 1 : numel(obj.names)
                nm  = obj.names{k};
                sz  = obj.sizes(k, :);
                lin = obj.offsets(k) : obj.offsets(k) + prod(sz) - 1;
                blk = reshape(w(lin), sz);
                if ~isempty(scaling)
                    blk = blk .* scaling.get(nm, sz(2));
                end
                parts.(nm) = blk;
            end
        end

        function idx = linear_indices(obj, name)
            % Column-major linear indices of the block inside w, as an
            % n_steps×n_cols matrix (same shape as the block).
            k  = obj.index_of(name);
            sz = obj.sizes(k, :);
            base = obj.offsets(k);
            idx = reshape(base : base + prod(sz) - 1, sz);
        end
    end

    methods (Access = private)
        function k = index_of(obj, name)
            k = find(strcmp(obj.names, char(name)), 1);
            assert(~isempty(k), 'DecisionLayout:unknown', ...
                'Unknown block "%s".', name);
        end
    end
end
