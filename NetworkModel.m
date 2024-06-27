classdef NetworkModel

    properties
        network;      % A BaseGraph object
        model;        % fluid model
        schedule;     % simulation schedule (wells and time steps)
        state0;       % initial state
        params = [];  % tunable parameters
        paramConfig;  % parameter configuration
        name = '';
    end

    methods
        function obj = NetworkModel(model, schedule, state0, network)
            % Network model for reservoir simulation.
            %
            % SYNOPSIS:
            %   networkmodel = NetworkModel(model, schedule, state0, network)
            %   networkmodel = NetworkModel(model, schedule, state0)
            %
            % REQIRED PARAMETERS:
            %   model    - Fluid model.
            %   schedule - Simulation schedule incl. wells and time steps.
            %   state0   - Initial state.
            %
            % OPTIONAL PARAMETERS:
            %   network    - Pre-made network. If not provided, a graph
            %                based on neighboring cells is used. 
            
            if nargin < 4
                % Convert grid direcly to network.
                obj.network = BaseGraph(model.operators.N);
                obj.model = model;
                coordNames = obj.coordNames();
                for i = 1:numel(coordNames)
                    obj.network.G.Nodes.(coordNames(i)) = model.G.cells.centroids(:,i);
                end

            else
                % Use given network.
                obj.network = BaseGraph(model.operators.N);
                obj.network.G = network;
                model.G.cells.centroids = obj.getNodeCoordinates();
            end

            % Model
            model = model.removeStateFunctionGroupings();
            model.G.cells.num = obj.network.numNodes;


            pv = mean(model.operators.pv)*ones(obj.network.numNodes, 1);
            T  = mean(model.operators.T)*ones(obj.network.numEdges, 1);
            model = model.setupOperators(model.G, ...
                        model.rock, ...
                        'trans', T,  ...
                        'neighbors', obj.network.G.Edges.EndNodes, ...
                        'porv', pv);

            obj.model = model;

            if isfield(model.rock, 'krscale')
                obj.model = obj.reimposeRelpermScaling();
            end

            obj.schedule = schedule;
            obj.state0 = state0;

            if ~isfield(obj.network.G.Nodes, 'x')
                obj = obj.assignCentroidCoordsToNodes();
            end
        end

        %------------------------------------------------------------------
        function model = updateModel(obj)
            model = obj.model.removeStateFunctionGroupings();
            model.G.cells.num = obj.network.numNodes;
            model.G.cells.centroids = obj.getNodeCoordinates();
            model = model.setupOperators(model.G, ...
                        model.rock, ...
                        'trans', model.operators.T,  ...
                        'neighbors', obj.network.G.Edges.EndNodes, ...
                        'porv', model.operators.pv);
        end

        %------------------------------------------------------------------
        function obj = rescalePoreVolumes(obj, origModel)
            % Rescale pore volumes so the sum matches the original model.
            origSumPV = sum(origModel.operators.pv);
            newPV = obj.model.operators.pv;
            obj.model.operators.pv = newPV * origSumPV/ sum(newPV);   
        end

        %------------------------------------------------------------------
        function state0 = expandState0(obj)
            % Increase size of state0 to match dimensions of graph.
            % Fill in with mean values. 
            newSize = struct('pressure', numNodes(obj.network), ...
                             'flux',     numEdges(obj.network), ...
                             's',        numNodes(obj.network));
            oldSize = struct('pressure', size(obj.state0.pressure,1), ...
                             'flux',     size(obj.state0.flux, 1), ...
                             's',        size(obj.state0.s, 1));
            fn = fieldnames(obj.state0);
            for i = 1:numel(fn)
                state0.(fn{i}) = [obj.state0.(fn{i}) ; 
                                  repmat(mean(obj.state0.(fn{i})(1:end,:)), ...
                                  newSize.(fn{i}) - oldSize.(fn{i}), 1)];
            end
        end

        %------------------------------------------------------------------
        function problem = getPackedSimulationProblem(obj, name)
            if nargin < 2
                name = obj.name;
            end
            problem = packSimulationProblem(obj.state0, ...
                                            obj.model, ...
                                            obj.schedule, ...
                                            name);
        end

        %------------------------------------------------------------------
        function sz = size(obj)
            % Returns number of nodes, edges and parameters
            sz = struct('Nodes', obj.network.numNodes, ...
                        'Edges', obj.network.numEdges);         
            if ~isempty(obj.params)
                sz.('Params') = obj.numParams();
            end
        end

        %=================================================================%
        %                        WELL UTILS                               %
        %=================================================================%
        
        function nW = numWells(obj)
            nW = numel(obj.schedule.control(1).W);
        end
        
        %------------------------------------------------------------------
        function wn = wellNodes(obj)
            wn = vertcat(obj.schedule.control(1).W.cells);
        end

        %=================================================================%
        %                          PARAMETERS                             %
        %=================================================================%
        
        function N = numParams(obj)
            N = cellfun(@(p)p.nParam, obj.params);
            N = sum(N);
        end

        %------------------------------------------------------------------
        function val = getParameterValue(obj, param)
            val = param.getParameterValue(obj);
        end
        
        %------------------------------------------------------------------
        function obj = setParameterValue(obj, param, val)
            obj = param.setParameterValue(obj, val);
        end

        %------------------------------------------------------------------
        function obj = setSingleParameterValue(obj, param, ix, val)
            p = obj.getParameterValue(param);
            p(ix) = val;
            obj = obj.setParameterValue(p);
        end
        
        %------------------------------------------------------------------
        function sens = getParamSensitivities(obj, params, ix)
            sens = cell2mat(cellfun(@(p) p.sens(ix), params, 'UniformOutput', false));
        end

        %------------------------------------------------------------------
        function np = nodeParams(obj)
            np = cell(0);
            for i = 1:numel(obj.params)
                if strcmp(obj.params{i}.mapTo, 'node')
                    np{end+1} = obj.params{i};
                end
            end
        end

        %------------------------------------------------------------------
        function np = nodeParamsExclState0(obj)
            np = cell(0);
            for i = 1:numel(obj.params)
                if strcmp(obj.params{i}.mapTo, 'node') && ...
                        ~strcmp(obj.params{i}.belongsTo, 'state0')
                    np{end+1} = obj.params{i};
                end
            end
        end

        %------------------------------------------------------------------
        function ep = edgeParams(obj)
            ep = cell(0);
            for i = 1:numel(obj.params)
                if strcmp(obj.params{i}.mapTo, 'edge')
                    ep{end+1} = obj.params{i};
                end
            end
        end

        %------------------------------------------------------------------
        function wp = wellParams(obj)
            wp = cell(0);
            for i = 1:numel(obj.params)
                if strcmp(obj.params{i}.mapTo, 'well')
                    wp{end+1} = obj.params{i};
                end
            end
        end

        %------------------------------------------------------------------
        function sp = state0Params(obj)
            sp = cell(0);
            for i = 1:numel(obj.params)
                if strcmp(obj.params{i}.belongsTo, 'state0')
                    sp{end+1} = obj.params{i};
                end
            end
        end

        %------------------------------------------------------------------
        function obj = fillinNodeParameters(obj, newnodes, varargin)
            % Fill in missing node parameters.
            % Iterate through, filling in parameters according to method
            % specified in param.fillStrategy.

            opt = struct('value', 'calcNow');
            opt = merge_options(opt, varargin{:});

            nodeparams = obj.nodeParams();
            for i = 1:numel(nodeparams)
                if strcmp(nodeparams{i}.belongsTo, 'state0') % does not work to set sw here
                    continue
                end
                val = obj.getParameterValue(nodeparams{i});
                for n = 1:numel(newnodes)
                    node = newnodes(n);
                    switch opt.value
                        case 'nan'
                            if node >= size(val)
                                val(node) = nan;
                            else
                                val = [val(1:node); nan; val(node+1:end)];
                            end
                        case 'remove'
                            val(node) = [];
                        case 'calcNow'  
                            val(node) = calcNodeParameter(obj.network, newnodes(n), val, nodeparams{i}.fillStrategy);
                    end
                end
                obj = obj.setParameterValue(nodeparams{i}, val);
            end
            % fixes nParam and calculates new boxLims
            obj.params = setupParameters(obj, obj.paramConfig);
        end

        %------------------------------------------------------------------
        function obj = fillinNanParameters(obj, rescaleNodes)
            % Find all nodes and edges with missing (nan) parameters.
            % Iterate through, filling in parameters according to methods
            % specified in param.fillStrategy.
            
            %--------------------------------------------------------------
            % Nodes
            nodeIDs = 1:obj.network.numNodes();

            % semi-local rescaling factor for pore volume
            pv0 = sum(obj.model.operators.pv(rescaleNodes), 'omitnan');

            nodeParameters = cell2mat(applyFunction(@(p) ...
                obj.getParameterValue(p), obj.nodeParamsExclState0));
            newNodes = nodeIDs(any(ismissing(nodeParameters),2));

            obj = obj.fillinNodeParameters(newNodes);

            % rescale pore volumes for new nodes and their neighbors
            obj.model.operators.pv(rescaleNodes) = obj.model.operators.pv(rescaleNodes) * pv0 / ...
                sum(obj.model.operators.pv(rescaleNodes));
            %--------------------------------------------------------------
            % Edges
            edgeIDs = 1:obj.network.numEdges();
            
            edgeParameters = cell2mat(applyFunction(@(p) ...
                obj.getParameterValue(p), obj.edgeParams));

            newEdges = edgeIDs(ismissing(edgeParameters));

            obj = obj.fillinEdgeParameters(newEdges);

        end

        %------------------------------------------------------------------
        function obj = fillinEdgeParameters(obj, newedges, varargin)
            % Fill in missing edge parameters.
            % Iterate through, filling in parameters according to method
            % specified in param.fillStrategy.

            opt = struct('value', 'calcNow');
            opt = merge_options(opt, varargin{:});

            edgeparams = obj.edgeParams();
            for i = 1:numel(edgeparams)
                val = obj.getParameterValue(edgeparams{i});
                for e = 1:numel(newedges)
                    edge = newedges(e);
                    switch opt.value
                        case 'nan'
                            if edge > size(val)
                                val(edge) = nan;
                            else
                                val = [val(1:edge); nan; val(edge+1:end)];
                            end
                        case 'remove'
                            val(edge) = [];
                        case 'calcNow'
                            val(edge) = calcEdgeParameter(obj.network, edge, ...
                                    val, edgeparams{i}.fillStrategy);  
                    end
                end
                obj = obj.setParameterValue(edgeparams{i}, val);
            end
        end            
       
        %=================================================================%
        %                           RELPERM                               %
        %=================================================================%

        function bool = hasRelperm(obj)
            % Check if model has relperm scaling
            bool = isfield(obj.model.rock, 'krscale');
        end

        %------------------------------------------------------------------
        function model = reimposeRelpermScaling(obj)
            % Reimpose relperm to match new dimension
            pts = obj.model.fluid.krPts;
            model = imposeRelpermScaling(obj.model, ...
                'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
                'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
                'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
        end

        %------------------------------------------------------------------
        function obj = expandRelperms(obj, newSize)
            oldSize = size(obj.model.rock.krscale.drainage.w, 1);
            obj.model.rock.krscale.drainage.w(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.drainage.ow(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.drainage.og(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.drainage.g(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.imbibition.w(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.imbibition.ow(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.imbibition.og(end+1:newSize,:) = nan(newSize-oldSize, 4);
            obj.model.rock.krscale.imbibition.g(end+1:newSize,:) = nan(newSize-oldSize, 4);
        end
        
        %=================================================================%
        %                            GRIDDING                             %
        %=================================================================%

        function d = dim(obj)
            d = size(obj.model.G.cells.centroids, 2);
        end

        %------------------------------------------------------------------
        function G = makeGrid(obj)
            % Make PEBI grid from graph
            % Currently assumes 2D
            coords = obj.getNodeCoordinates();
            if obj.dim == 3
                x = coords(:,1);
                y = coords(:,2);
                z = coords(:,3);
                assert(all(z == z(1)));
            else
                x = coords(:,1);
                y = coords(:,2);
            end
            X   = [x, y];
            bix = convhull(X);
            G   = clippedPebi2D(X, X(bix,:));
        end

        %=================================================================%
        %                          PLOTTING                               %
        %=================================================================%
        
        function plot(obj, varargin)
            plotNetworkModel(obj, varargin{:});
        end

        %------------------------------------------------------------------
        function plotCircular(obj, varargin)
            plotNetworkModelCircular(obj, varargin{:});
        end

        %=================================================================%
        %                        COORDINATES                              %
        %=================================================================%
            
        function varargout = getCentroidCoordinates(obj)
            % [x, y] = obj.getCentroidCoordinates() for 2D
            % [x, y, z] = obj.getCentroidCoordinates() for 3D
            nOutputs = nargout;
            varargout = cell(1,nOutputs);
            for i = 1:nOutputs
                varargout{i} = obj.model.G.cells.centroids(:, i);
            end
        end

        %------------------------------------------------------------------
        function coords = getNodeCoordinates(obj)
            coords(:,1) = obj.network.G.Nodes.x;
            coords(:,2) = obj.network.G.Nodes.y;
            if any(strcmp(obj.network.G.Nodes.Properties.VariableNames,'z'))
                coords(:,3) = obj.network.G.Nodes.z;
            end
        end

        %------------------------------------------------------------------
        function cn = coordNames(obj)
           dim = size(obj.model.G.cells.centroids, 2);
           if dim == 2
               cn = ["x", "y"];
           else
               cn = ["x", "y", "z"];
           end
        end

        %------------------------------------------------------------------
        function obj = assignCentroidCoordsToNodes(obj)
            cn = obj.coordNames();
            if numel(cn) == 2
                [x,y] = obj.getCentroidCoordinates();
                obj.network.G.Nodes.('x') = x;
                obj.network.G.Nodes.('y') = y;
            elseif numel(cn) == 3
                [x,y,z] = obj.getCentroidCoordinates();
                obj.network.G.Nodes.('x') = x;
                obj.network.G.Nodes.('y') = y;
                obj.network.G.Nodes.('z') = z;
            end
        end

        %------------------------------------------------------------------
        function x = x(obj)
            %x = obj.model.G.cells.centroids(:,1);
            x = obj.network.G.Nodes.x;
        end

        function y = y(obj)
            %y = obj.model.G.cells.centroids(:,2);
            y = obj.network.G.Nodes.y;
        end

        function z = z(obj)
            %z = obj.model.G.cells.centroids(:,3);
            z = obj.network.G.Nodes.z;
        end

        %=================================================================%
        %                      GRAPH MODIFICATION                         %
        %=================================================================%

        function [obj, nodeID] = addNode(obj)
            [obj.network, nodeID] = obj.network.addNode();
            if obj.hasRelperm
                obj = obj.expandRelperms(obj.network.numNodes);
            end
            obj = obj.fillinNodeParameters(nodeID, 'value', 'nan');
        end

        %------------------------------------------------------------------
        function obj = removeNode(obj, nodeID)
            obj.network = obj.network.removeNode(nodeID);
            obj = obj.fillinNodeParameters(nodeID, 'value', 'remove');
        end

        %------------------------------------------------------------------
        function [obj, edgeID] = addEdge(obj, s, t)
            [obj.network, edgeID] = obj.network.addEdge(s,t);
            obj = obj.fillinEdgeParameters(edgeID, 'value', 'nan');
        end

        %------------------------------------------------------------------
        function obj = removeEdge(obj, s, t)
            [obj.network, edgeID] = obj.network.removeEdge(s,t);
            obj = obj.fillinEdgeParameters(edgeID, 'value', 'remove');
        end

        %=================================================================%
        %                           CHECKS                                %
        %=================================================================%

        function checkDimensions(obj)
            assert(obj.model.G.cells.num == obj.network.numNodes, ...
                    'Number of cells and nodes not equal.')
            assert(size(obj.model.operators.N, 1) == obj.network.numEdges, ...
                    'Number of neighbor cells and edges not equal.')
        end
    end
end
