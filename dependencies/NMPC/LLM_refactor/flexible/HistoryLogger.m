classdef HistoryLogger < handle
    % HistoryLogger  Append-only array of solve records.

    properties (SetAccess = private)
        records = struct([])
    end

    methods
        function append(obj, rec)
            if isempty(obj.records)
                obj.records = rec;
            else
                obj.records(end + 1) = rec;
            end
        end

        function n = count(obj)
            n = numel(obj.records);
        end

        function clear(obj)
            obj.records = struct([]);
        end
    end
end
