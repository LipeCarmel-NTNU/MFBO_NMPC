classdef FminconBackend < handle
    % FminconBackend  Wraps fmincon behind a single problem-struct interface.

    properties (SetAccess = private)
        options
    end

    methods
        function obj = FminconBackend(args)
            arguments
                args.options = []
            end
            if isempty(args.options)
                obj.options = optimoptions('fmincon', ...
                    'Display','off','Algorithm','sqp', ...
                    'MaxFunEvals',Inf,'MaxIterations',1000, ...
                    'StepTolerance',1e-9,'OptimalityTolerance',1e-6, ...
                    'ScaleProblem',true);
            else
                obj.options = args.options;
            end
        end

        function [w_opt, fval, flag] = solve(obj, problem)
            opts = problem.options;
            if isempty(opts)
                opts = obj.options;
            end
            [w_opt, fval, flag] = fmincon( ...
                problem.objfun, problem.w0, ...
                problem.A,   problem.b, ...
                problem.Aeq, problem.beq, ...
                problem.wL,  problem.wU, ...
                problem.nonlcon, opts);
        end
    end
end
