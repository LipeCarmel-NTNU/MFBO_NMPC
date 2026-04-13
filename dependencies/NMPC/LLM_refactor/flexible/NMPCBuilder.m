classdef NMPCBuilder < handle
    % NMPCBuilder  Fluent assembler that produces a validated NMPC.
    %
    %   nmpc = NMPCBuilder() ...
    %       .withModel(ContinuousModel(@f, nx, nu, RK4(Ts))) ...
    %       .withHorizon(p, m) ...
    %       .withScaling('x', x_scale).withScaling('u', u_scale) ...
    %       .addConstraint(StateBounds(Ymin, Ymax)) ...
    %       .addConstraint(InputBounds(umin, umax)) ...
    %       .addCost(QuadraticTracking(Q, R, y_sp, u_sp)) ...
    %       .enableLog() ...
    %       .build();

    properties (Access = private)
        model_      = []
        p_          = []
        m_          = []
        ny_         = []
        scaling_
        constraints_ = {}
        costs_       = {}
        backend_     = []
        log_enabled_ = false
    end

    methods
        function obj = NMPCBuilder()
            obj.scaling_ = Scaling();
        end

        function obj = withModel(obj, model)
            obj.model_ = model;
        end

        function obj = withHorizon(obj, p, m)
            obj.p_ = p;
            obj.m_ = m;
        end

        function obj = withOutputs(obj, ny)
            obj.ny_ = ny;
        end

        function obj = withScaling(obj, name, factor)
            obj.scaling_.set(name, factor);
        end

        function obj = addConstraint(obj, c)
            obj.constraints_{end + 1} = c;
        end

        function obj = addCost(obj, c)
            obj.costs_{end + 1} = c;
        end

        function obj = withBackend(obj, backend)
            obj.backend_ = backend;
        end

        function obj = enableLog(obj, tf)
            if nargin < 2
                tf = true;
            end
            obj.log_enabled_ = logical(tf);
        end

        function nmpc = build(obj)
            % Validate.
            assert(~isempty(obj.model_),  'NMPCBuilder:model',  'Missing model.');
            assert(~isempty(obj.p_),      'NMPCBuilder:p',      'Missing horizon p.');
            assert(~isempty(obj.m_),      'NMPCBuilder:m',      'Missing horizon m.');
            assert(obj.m_ <= obj.p_,      'NMPCBuilder:horizon','Require m <= p.');
            assert(~isempty(obj.costs_),  'NMPCBuilder:cost',   'At least one cost required.');

            nx = obj.model_.nx;
            nu = obj.model_.nu;
            ny = obj.ny_;
            if isempty(ny); ny = nx; end

            % Build layout. x and u are always present.
            layout = DecisionLayout();
            layout.add('x', obj.p_ + 1, nx);
            layout.add('u', obj.m_,     nu);

            % Let constraints register extra blocks (slacks, aux).
            for k = 1 : numel(obj.constraints_)
                c = obj.constraints_{k};
                if ismethod(c, 'register')
                    c.register(layout);
                end
            end
            % And costs, in case a cost wants to contribute blocks.
            for k = 1 : numel(obj.costs_)
                c = obj.costs_{k};
                if ismethod(c, 'register')
                    c.register(layout);
                end
            end
            layout.freeze();

            % Continuity is always required for multiple shooting.
            constraints = [obj.constraints_, {Continuity()}];

            backend = obj.backend_;
            if isempty(backend)
                backend = FminconBackend();
            end

            logger = [];
            if obj.log_enabled_
                logger = HistoryLogger();
            end

            cfg = struct( ...
                'model',       obj.model_, ...
                'p',           obj.p_, ...
                'm',           obj.m_, ...
                'ny',          ny, ...
                'layout',      layout, ...
                'scaling',     obj.scaling_, ...
                'constraints', {constraints}, ...
                'costs',       {obj.costs_}, ...
                'backend',     backend, ...
                'log_enabled', obj.log_enabled_, ...
                'logger',      logger);

            nmpc = NMPC(cfg);
        end
    end
end
