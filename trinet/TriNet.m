classdef TriNet < NetworkModel
    % TriNet - Triangulation-based network model.
    %
    % SYNOPSIS:
    %   trinet = TriNet(model, schedule, state, 'pn1', pv1, ...)
    %   
    % DESCRIPTION:
    %   Constructs a network model based on a triangulation of well nodes
    %   and selected points along the boundary. 
    %   The class also holds methods for refining the graph.
    %
    % REQUIRED PARAMETERS:
    %   model    - Fluid model.
    %   schedule - Simulation schedule incl. wells and time steps.
    %   state0   - Initial state.
    %
    % OPTIONAL PARAMETERS:
    %   Any optional parameters in delaunayGraph/delaunayTriangulate.

    properties
        tri
    end
    
    methods

        function obj = TriNet(model, schedule, state0, varargin)
            [tri, G, wellNodes, state0] = delaunayGraph(model, schedule, state0, varargin{:});
            obj = obj@NetworkModel(model, schedule, state0, G);
            obj.tri = tri;

            wellNodes = num2cell(wellNodes);
            for i = 1:numel(obj.schedule.control)
                [obj.schedule.control(i).W.cells] = wellNodes{:};
            end

            obj = obj.rescalePoreVolumes(model);
        end

        %==================================================================
        %                    TRIANGULATION UTILITIES                      %
        %==================================================================
        
        function n = numTri(obj)
            % Number of triangles in the triangulation.
            n = size(obj.tri, 1);
        end
    
        %------------------------------------------------------------------
        function nodes = getTriangleNodes(obj, triID)
            % Return the three corner nodes of given triangle. 
            n1 = obj.tri(triID, 1);
            n2 = obj.tri(triID, 2);
            n3 = obj.tri(triID, 3);
            nodes = [n1, n2, n3];
        end

        %------------------------------------------------------------------
        function edges = getTriangleEdges(obj,triID)
            % Return the three edges/sides of given triangle.
            nodes = obj.getTriangleNodes(triID);
            e1 = findedge(obj.network.G, nodes(1), nodes(2));
            e2 = findedge(obj.network.G, nodes(2), nodes(3));
            e3 = findedge(obj.network.G, nodes(1), nodes(3));
            edges = [e1, e2, e3];
        end
        
        %------------------------------------------------------------------
        function triangles = getSurroundingTriangles(obj, nodeID)
            % Find the triangles around a given node / triangles
            % where the given node is a corner.
            triangles = find(any(obj.tri == nodeID, 2));
        end

        %------------------------------------------------------------------
        function triangles = getEdgeTriangles(obj, s, t)
            % Find the triangles that a given edge (s,t) is part of.  
            % (both end nodes are in the triangle)
            triangles = find(any(ismember(obj.tri, s), 2) & ...
                             any(ismember(obj.tri, t), 2));
        end

        %------------------------------------------------------------------
        function nodes = findEdgeOppositeNodes(obj, s, t)
            % Find the third and last nodes of the triangles
            % an edge is part of.
            triangles = obj.getEdgeTriangles(s,t);
            triangleNodes = obj.tri(triangles,:);
            nodes = triangleNodes(~ismember(triangleNodes, [s,t]));
        end

        %------------------------------------------------------------------
        function neighborTriangle = triangleNeighborOnEdge(obj, triangle, s, t)
            % Find the triangle that shares the edge with given triangle 

            if nargin < 4 % Edge ID given instead of s,t
                [s,t] = obj.getEndNodes(s);
            end
            trianglesSharingEdge = find(any(ismember(obj.tri, s), 2) & ...
                                        any(ismember(obj.tri, t), 2));
            % Choose the one not equal to known triangle.
            neighborTriangle = trianglesSharingEdge(trianglesSharingEdge ~= triangle);
        end
        
        %==================================================================
        %                  GRAPH MODIFICATION ROUTINES                    %
        %==================================================================

        function [obj, rescaleNodes] = refineGraph(obj, opt)
            % Main graph refinement routine.
            rescaleNodes = 1:obj.network.numNodes; % fallback for edges and nodes case
            switch opt.refineFrom
                case 'triangles'
                    assert(~isempty(obj.tri), ...
                        'Cannot refine from triangles since triangulation is not defined');
                    selected = obj.selectTriangles(opt);
                    [obj, rescaleNodes] = refineTriangles(obj, selected);

                case 'edges'
                    selected = obj.selectEdges(opt);
                    edges = cell(numel(selected),1);
                    for i = 1:numel(selected)
                        [s,t] = obj.network.getEndNodes(selected(i));
                        edges{i} = [s,t];
                    end
                    for i = 1:numel(selected)
                        obj = obj.refineEdge(edges{i});
                    end
                case 'nodes'
                    selected = obj.selectNodes(opt);
                    obj = obj.refineNodes(selected);
                otherwise
                    error('Unsupported graph refinement method: %s', opt.refineFrom)
            end
        end

        %------------------------------------------------------------------
        function obj = refineAllTriangles(obj)
            % Uniform refinement of triangulation graph
            numNodesBefore = obj.network.numNodes;
            obj = refineTriangles(obj, 1:numTri(obj)); % Refine all triangles

            obj = obj.fillinEdgeParameters(1:obj.network.numEdges);
            obj = obj.fillinNodeParameters(numNodesBefore+1:obj.network.numNodes);
            
            obj.model = obj.updateModel();
            obj.state0 = obj.expandState0();
        end

        %------------------------------------------------------------------
        %                       Node refinement                           %
        %------------------------------------------------------------------
        function nodeSelect = selectNodes(obj, opt)
            % Select nodes to refine the graph around.
            %
            % SYNOPSIS:
            %   nodeSelect = obj.selectNodes(opt)
            % 
            % DESCRIPTION:
            %   Based on sensitivities, select which nodes needs refinement.
            %   Set 'selectWith' to 'anyAbove' to select the nodes with at 
            %   least one sensitivity above maxSens, or 'fraction' to
            %   select the top 'fraction' percent. 
            %   With 'anyAbove', if more nodes than maxNodes are selected, 
            %   a random subset is chosen.
            %
            % PARAMETERS:
            %   opt   - struct with fields:
            %       'selectWith' - How to select nodes. Supported options
            %                      are 'anyAbove' and 'fraction'.
            %       'maxSens'    - Threshold for the sensitivities
            %       'maxFrac'    - Maximum fraction of nodes in selection
            %        
            % RETURNS:
            %   Array of selected nodeIDs.
            
            nodeIDs = 1:obj.network.numNodes();

            S = obj.getParamSensitivities(obj.nodeParams, nodeIDs);

            maxNodes = round(opt.maxFrac * obj.network.numNodes());

            switch opt.selectWith
                case 'anyAbove'
                    nodeSelect = nodeIDs(S > opt.maxSens);
                    if numel(nodeSelect) > maxNodes
                        subset = randi([1, numel(nodeSelect), [1, maxNodes]]); 
                        nodeSelect = nodeSelect(subset);  
                    end
                case 'fraction'
                    [largest_vals, nodeSelect] = maxk(S, maxNodes);
                otherwise
                    error('Unsupported node selection method: %s', opt.selectWith)
            end
        end

        %------------------------------------------------------------------
        function obj = refineNodes(obj, nodeIDs)
            % Graph refinement from node. 
            % Refines surrounding triangles.

            % Find surrounding triangles.
            triSelect = [];
            for i = 1:numel(nodeIDs)
                node = nodeIDs(i);
                surroundTri = obj.getSurroundingTriangles(node);
                triSelect = [triSelect; surroundTri];
            end
            triSelect = unique(triSelect);

            % Refine surrounding triangles.
            obj = refineTriangles(obj, triSelect);
        end

        %------------------------------------------------------------------
        %                        Edge refinement                          %
        %------------------------------------------------------------------
        function edgeSelect = selectEdges(obj, opt)
            % Select edges to refine the graph around from sensitivities.
            % Selects edges with at least one sensitivity above maxSens.
            % If more edges than maxEdges are selected, a random
            % subset is chosen.
            %
            % PARAMETERS:
            %   opt - Struct with fields:
            %       'selectWith' - How to select edges. Supported options
            %                      are 'anyAbove' and 'fraction'.
            %       'maxSens'    - Threshold for the sensitivities
            %       'maxFrac'    - Maximum fraction of edges in selection
            % 
            % RETURNS:
            %   Array of selected edgeIDs.

            edgeIDs = 1:obj.network.numEdges();

            S = obj.getParamSensitivities(obj.edgeParams, edgeIDs);
            maxEdges = round(opt.maxFrac * obj.network.numEdges());

            switch opt.selectWith
                case 'anyAbove'
                    edgeSelect = edgeIDs(S > opt.maxSens);
                    if numel(edgeSelect) > maxNodes
                        subset = randi([1, numel(edgeSelect), [1, maxNodes]]);
                        edgeSelect = nodeSelect(subset);
                    end
                case 'fraction'
                    [largest_vals, edgeSelect] = maxk(S, maxEdges);
                otherwise
                    error('Unsupported edge selection method: %s', opt.selectWith)
            end
        end

        %------------------------------------------------------------------
        function [obj, nodeID] = splitEdge(obj, s, t, nodeID)
            % Split edge in two, adding node in middle.
            %   o---o  ->  o-x-o
            if nargin < 4
                [obj, nodeID] = obj.addNode();
            end
            obj.network.G.Nodes(nodeID, obj.coordNames).Variables = mean( ...
                obj.network.G.Nodes([s, t], obj.coordNames).Variables);
            obj = obj.addEdge(s, nodeID);
            obj = obj.addEdge(t, nodeID);
            obj = obj.removeEdge(s,t);
        end

        %------------------------------------------------------------------
        function obj = refineEdge(obj, edge)
            % Graph refinement from edge.
            % Splits edge by adding node in middle.
            % Connect new node to the third nodes of the triangles
            % the old edge is in. 
            %   o          o 
            %  / \   –>   /|\
            % o---o      o-x-o
            [s, t] = deal(edge(1), edge(2));
            nodes = obj.findEdgeOppositeNodes(s,t);
            [obj, nodeID] = obj.splitEdge(s,t);
            for i = 1:numel(nodes)
                obj = obj.addEdge(nodeID, nodes(i));
                % Update tri
                oldTri = find(any(ismember(obj.tri, s), 2) & ...
                              any(ismember(obj.tri, t), 2) & ...
                              any(ismember(obj.tri, nodes(i)), 2));
                obj.tri(oldTri, :) = [];
                obj.tri(end+1, :) = [s, nodeID, nodes(i)];
                obj.tri(end+1, :) = [t, nodeID, nodes(i)];
            end
        end

        %------------------------------------------------------------------
        %                        Triangle refinement                      %
        %------------------------------------------------------------------
        function triSelect = selectTriangles(obj, opt)
            % Select triangles for graph refinement.
            %
            % SYNOPSIS:
            %   triSelect = g.selectTriangles()
            %       uses default choice (any sensitivity above thres)
            %   triSelect = g.selectTriangles('selectWith', 'fraction', 'maxFrac', 0.2)
            %       selects top 20% with largest sensitivity.
            %
            % DESCRIPTION: 
            %   Select triangles for graph refinement based on
            %   sensitivities of the triangle's nodes and edges. 
            %   Assuming we have sensitivities mapped to the nodes and
            %   edges of the graph, we map them to triangles, and choose
            %   triangles to refine based on the mapped values. 
            %
            % PARAMETERS:
            %   opt   - struct with fields:
            %       'sensMap'       - Mapping of sensitivities from nodes and edges
            %                         to triangles. 
            %                         Options: 'max', 'sum', 'mean'.
            %       'selectWith'    - How to select the triangles. Supported 
            %                         options are 'anyAbove' and 'fraction'.
            %       'maxSens'       - Threshold for the sensitivities.
            %       'maxFrac'       - Maximum fraction of triangles allowed
            %                         allowed in selection.
            %                       
            % RETURNS:
            %   triSelect       - Array of selected triangles (indices).

            switch opt.sensMap
                case 'max'
                    % Maximum of all sensitivities
                    indMap = @(triID) max( ...
                        max(obj.getParamSensitivities(obj.edgeParams, [obj.getTriangleEdges(triID)]), [], 'all'), ...
                        max(obj.getParamSensitivities(obj.nodeParams, [obj.getTriangleNodes(triID)]), [], 'all'));
                case 'sum'
                    % Sum of all sensitivities
                    indMap = @(triID) sum([obj.getParamSensitivities(obj.edgeParams, [obj.getTriangleEdges(triID)]), ...
                                            obj.getParamSensitivities(obj.nodeParams, [obj.getTriangleNodes(triID)])], 'all');
                case 'mean'
                    % 
                    indMap = @(triID) mean([obj.getParamSensitivities(obj.edgeParams, [obj.getTriangleEdges(triID)]), ...
                                             obj.getParamSensitivities(obj.nodeParams, [obj.getTriangleNodes(triID)])], 'all');
                otherwise
                    error('Unsupported triangle sensitivity mapping: %s', opt.sensMap)
            end
            
            triIDs = 1:numTri(obj);
            val = arrayfun(@(triID) indMap(triID), triIDs');  % refinement indicator values for all triangles
            maxTri = round(opt.maxFrac * numTri(obj));        % number of triangles to select
            
            switch opt.selectWith
                case 'anyAbove'
                    % Select triangles with any sensitivity above threshold
                    % If more triangles than allowed are selected, we draw
                    % a random subset.
                    triSelect = triIDs(val > opt.maxSens);
                    if numel(triSelect) > maxTri
                        subset = randi([1, numel(triSelect)], [1, maxTri]); % maxTri random triangle indices.
                        triSelect = triSelect(subset);
                    end
                case 'fraction'
                    % Select the triangles with largest sensitivity value.
                    [largest_vals, triSelect] = maxk(val, maxTri);     
                otherwise
                    error('Unsupported triangle selection method: %s', opt.selectWith)
            end
        end
    end
end
