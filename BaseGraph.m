 classdef BaseGraph

    % Base class implementing basic graph functionality.
    %
    % The class is abstract and BaseGraph objects can not be
    % constructed directly.
    %
    % DESCRIPTION:
    %   This class serves as a base for more specialized simulation
    %   graphs and aims to implement utilities that are common for 
    %   all types of graphs. 
    %   The property G is an instance of a MATLAB graph, and many 
    %   class methods passes to methods of the built-in MATLAB 
    %   graph class.

    properties
        G;  % A MATLAB graph object
    end

    methods
        function obj = BaseGraph(N)
            obj.G = graph(N(:,1), N(:,2));
        end

        %------------------------------------------------------------------
        function n = numNodes(obj)
            % Number of nodes in the graph.
            n = numnodes(obj.G);
        end

        %------------------------------------------------------------------
        function n = numEdges(obj)
            % Number of edges in the graph.
            n = numedges(obj.G);
        end
        
        %------------------------------------------------------------------
        function ns = getNeighbors(obj, nodeID)
            % Return nodeIDs for the neighbor nodes of given node.
            ns = obj.G.neighbors(nodeID);
        end

        %------------------------------------------------------------------
        function [s,t] = getEndNodes(obj, edgeID)
            % Return nodeIDs for the two nodes connected by the given edge.
            s = obj.G.Edges.EndNodes(edgeID, 1);
            t = obj.G.Edges.EndNodes(edgeID, 2);
        end    
        
        %------------------------------------------------------------------
        function edges = getEdgesFromNode(obj, nodeID)
            % Get list of edgeIDs for edges from the given node.
            neighbors = obj.getNeighbors(nodeID);
            edges = zeros(numel(neighbors),1);
            for n = 1:numel(neighbors)
                % Find edgeID for edge to neighbor.
                edges(n) = findedge(obj.G, nodeID, neighbors(n));
            end
        end

        %==================================================================
        %----------- MATRIX REPRESENTATION --------------------------------
        %==================================================================
        
        function A = adjacency(obj)
            % Get the adjacency matrix of the graph.
            % Size: #Nodes x #Nodes
            %   A_ij = 1 if there is an edge between node i and node j,
            %          0 otherwise
            A = obj.G.adjacency;
        end

        %------------------------------------------------------------------
        function I = incidence(obj)
            % Get the incidence matrix of the graph.
            % Size: #Nodes x #Edges
            %   B_ij = 1 if node i and edge j are incident (node i is on edge j), 
            %          0 otherwise
            I = obj.G.incidence;
        end

        %------------------------------------------------------------------
        function L = laplacian(obj)
            % Get the Laplacian matrix of the graph.
            % L = I^T * I, I being the incidence matrix.
            L = obj.G.laplacian;
        end
        
        %------------------------------------------------------------------
        function C = sharesNeighbors(obj)
            % Returns list of node pairs who share at least one 
            % neighbor node.
            % 
            % EXAMPLE:
            %   If there are two edges (1,2) and (2,3), then 1 and 3 
            %   share the neighbor 2 and will appear as a row in C.
            %
            A = obj.adjacency;
            % A^T * A, since A is symmetric for undirected graphs, we
            % skip the transpose, and use A*A.
            % The upper-triangular part is extracted to avoid duplicates.
            C = triu(A*A); 
            % Extract indices of non-zero entries. (corr. to nodeIDs)
            [row, col] = find(C);
            C = [row, col];
            C = C(C(:,1)~=C(:,2),:); % remove diagonal elements.
        end

        %------------------------------------------------------------------
        function N = getCommonNeighbors(obj, n1, n2)
            % Returns a list of nodes having edges
            % to both n1 and n2, that is, the common 
            % neighbors of n1 and n2.
            %
            % EXAMPLE:
            %   If there are two edges (1,2) and (2,3), then
            %   obj.getCommonNeighbors(1,3) will return [2].
            %
            N1 = obj.getNeighbors(n1);
            N2 = obj.getNeighbors(n2);
            N = intersect(N1, N2);
        end
        
        %------------------------------------------------------------------
        %---------- GRIDDING ----------------------------------------------
        %------------------------------------------------------------------
        function G = makeGrid(obj)
            % Make PEBI grid from triangulation
            % Currently assumes 2D
            if any(strcmp('z', obj.G.Nodes.Properties.VariableNames))
                assert(all(obj.G.Nodes.z == obj.G.Nodes.z(1)));
            end
            x   = [obj.G.Nodes.x, obj.G.Nodes.y];
            bix = convhull(x);
            G   = clippedPebi2D(x, x(bix,:));
        end
        
        %------------------------------------------------------------------
        %----------- CHECKS -----------------------------------------------
        %------------------------------------------------------------------
        function removeNodeCheck(obj, nodeID)
            % Check that node exists.
            if ~findnode(obj.G, nodeID)
                warning('Attempted to remove non-existing node. Graph is unchanged.')
            end
        end

        %------------------------------------------------------------------
        function checkEdge(obj, s, t)
            % Use when adding new edge. 
            % Checks if end nodes exist, and
            % checks if edge already exists. (duplicates not allowed)
            assert(findnode(obj.G,s) && findnode(obj.G,t), 'Nodes do not exist.')
            assert(~findedge(obj.G,s,t), 'Edge already exists.')
        end

        %------------------------------------------------------------------
        %----------- ADD AND REMOVE NODES AND EDGES -----------------------
        %------------------------------------------------------------------
        
        function [obj, nodeID] = addNode(obj)
            obj.G = addnode(obj.G, 1);
            nodeID = obj.numNodes; % added last
        end

        %------------------------------------------------------------------
        function obj = removeNode(obj, nodeID)
            obj.removeNodeCheck(nodeID);
            obj.G = rmnode(obj.G, nodeID);
        end
        
        %------------------------------------------------------------------
        function [obj, edgeID] = addEdge(obj, s, t)
            obj.G = addedge(obj.G, s, t);
            edgeID = findedge(obj.G, s, t);
        end

        %------------------------------------------------------------------
        function [obj, edgeID] = removeEdge(obj, s, t)
            % Remove the edge from node s to node t, or 
            % issue a warning if the edge does not exist.
            edgeID = findedge(obj.G, s, t);
            if ~edgeID
                warning('Attempted to remove non-existing edge (%d, %d).', s, t)
            else
                obj.G = rmedge(obj.G, s, t);
            end
        end



    end
 end
