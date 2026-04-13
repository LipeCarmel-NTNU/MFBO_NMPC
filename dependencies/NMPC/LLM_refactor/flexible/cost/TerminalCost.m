classdef TerminalCost < handle
    % TerminalCost  ‖x_N - x_sp‖_P².

    properties (SetAccess = private)
        P
        x_sp
    end

    methods
        function obj = TerminalCost(P, x_sp)
            arguments
                P    double
                x_sp (1,:) double
            end
            obj.P    = P;
            obj.x_sp = x_sp;
        end

        function J = evaluate(obj, parts, ~)
            e = (parts.x(end, :) - obj.x_sp).';
            J = e.' * obj.P * e;
        end
    end
end
