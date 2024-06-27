function [trinet, rescaleNodes] = refineTriangles(trinet, selected)   
    % Refine the TriNet by splitting selected triangles.
    %
    % SYNOPSIS: trinet = refineTriangles(trinet, selected)
    %
    % DESCRIPTION:
    %   Graph refinement inspired by FEM red-green mesh refinement.
    %   The selected triangles are refined in the red step by split in 
    %   four. Then the hanging nodes are removed in the green step.
    %
    % REQUIRED PARAMETERS: 
    %   trinet    - TriNet (NetworkModel)
    %   selected  - List of selected triangles (row numbers in tri matrix).
    %
    % RETURNS:
    %   trinet - Refined TriNet.
    %   rescaleNodes - Nodes of selected triangles and their neighbors.

    % Initialize status struct. 
    status = struct('toRefine',  selected, ...  % Triangles to be refined.
                    'inspected', [], ...        % Inspected triangles (for recursive step).
                    'refined',   [], ...        % Refined triangles.
                    'flagged',   [], ...        % Flagged (red) triangles.
                    'splitNodes',  []); % 2D array holding rows [s, t, splitNodeID].
    
    if isempty(status.toRefine)
        warning('No triangles were selected for refinement.')
        rescaleNodes = [];
        return
    end
    
    %---RED STEP-----
    tri0 = trinet.tri;
    while ~isempty(status.toRefine)
        [trinet, status] = refineSingleTriangle(trinet, tri0, status.toRefine(1), status);
    end

    rescaleNodes = [status.splitNodes(:,3)];
    for i = 1:numel(selected)
        rescaleNodes = [rescaleNodes; tri0(i,:)'];
    end
    for i = 1:size(status.flagged,1)
        triangle = status.flagged(i, 1);
        rescaleNodes = [rescaleNodes; tri0(triangle,:)'];
    end
    rescaleNodes = unique(rescaleNodes);

    %---GREEN STEP----
    [trinet, status] = refineFlaggedTriangles(trinet, tri0, status);
end

function [trinet, status] = refineSingleTriangle(trinet, tri0, triID, status)
    % Refinement of a single triangle.
    %
    % SYNOPSIS: 
    %   [trinet, status] = refineSingleTriangle(trinet, tri0, triID, status)
    %
    % DESCRIPTION:
    %   The triangle at tri0(triID,:) is refined by splitting in four. 
    %   The status of the neighbouring triangles are checked and:
    %       - If neighbor triangle should be refined, the current triangle
    %         is marked as inspected and the neighbour is refined first in
    %         a recursive manner.
    %       - If neighbor triangle is not selected for refinement, it is
    %         flagged. We keep track of how many times N a triangle is
    %         flagged in status.flagged. If triangle i is flagged N times,
    %         status.flagged has a row [i, N].
    %
    % REQUIRED PARAMETERS:
    %   trinet  - TriNet (NetworkModel with triangulation)
    %   tri0    - Triangulation before any refinement (so selected rows match).
    %   triID   - The triangle to refine, find corners by tri0(triID,:).
    %   status  - Struct with info about refinements, flagging, edge splits 
    %             etc. (See refineTriangles)
    %
    % RETURNS:
    %   trinet  - Updated TriNet
    %   status  - Updated status

    triangle = tri0(triID,:);
    triNodes = sort(triangle); % ascending order
    triEdges = [triNodes(1), triNodes(2);
                triNodes(2), triNodes(3);
                triNodes(1), triNodes(3)];
    status.inspected = [status.inspected; triID];

    newnodes = zeros(1,3);
    for j = 1:3
        edge = triEdges(j,:);
        s = edge(1);
        t = edge(2);

        % Find triangle on opposite side of edge.
        trianglesSharingEdge = find(any(ismember(tri0, s), 2) & ...
                                    any(ismember(tri0, t), 2));
        oppTri = trianglesSharingEdge(trianglesSharingEdge ~= triID);


        if ~isempty(oppTri) && ismember(oppTri, status.inspected)
            % If inspected but not refined, we make the split node now. 
            if ~ismember(oppTri, status.refined)
                [trinet, newnode] = trinet.addNode();
                trinet = trinet.splitEdge(s, t, newnode);
            % If already refined, find and use the existing split node.
            else
                newnode = status.splitNodes(all(status.splitNodes(:,1:2)==edge,2), 3);
            end

        elseif ~isempty(oppTri) && ismember(oppTri, status.toRefine)
            % If neighbour triangle is to be refined, do it now and come
            % back to current triangle.
            [trinet, status] = refineSingleTriangle(trinet, tri0, oppTri, status);
            newnode = status.splitNodes(all(status.splitNodes(:,1:2)==edge,2), 3);
            
        else
            % Add splitting node
            [trinet, newnode] = trinet.addNode();
            
            % Split the edge.
            trinet = trinet.splitEdge(edge(1), edge(2), newnode);
        end
        
        % Flag neighbour triangle if it is not selected for refinement.
        if ~isempty(oppTri) && ...
           ~ismember(oppTri, status.refined) && ...
           ~ismember(oppTri, status.toRefine)
            status = flagTriangle(status, oppTri);
        end
        
        newnodes(j) = newnode; % Add to current triangle's new nodes
        status.splitNodes = [status.splitNodes; s, t, newnode]; % Keep split node for lookup later
    end

    % Add edges between the new nodes
    trinet = trinet.addEdge(newnodes(1), newnodes(2));
    trinet = trinet.addEdge(newnodes(2), newnodes(3));
    trinet = trinet.addEdge(newnodes(1), newnodes(3));

    % Update tri
    newRow = find(ismember(trinet.tri, triangle, 'rows'));
    trinet.tri(newRow,:) = []; % Remove old triangle
    trinet.tri(end+1,:) = [triNodes(1), newnodes(1), newnodes(3)]; % corner 1
    trinet.tri(end+1,:) = [triNodes(2), newnodes(1), newnodes(2)]; % corner 2
    trinet.tri(end+1,:) = [triNodes(3), newnodes(2), newnodes(3)]; % corner 3
    trinet.tri(end+1,:) = [newnodes(1), newnodes(2), newnodes(3)]; % middle
        
    % Move triangle from toRefine to refined.
    status.toRefine = status.toRefine(status.toRefine ~= triID);
    status.refined = [status.refined; triID];
end

function [trinet, status] = refineFlaggedTriangles(trinet, tri0, status)
    % Green step: Refine flagged triangles to remove hanging nodes.
    %   
    % SYNOPSIS:
    %   [trinet, status] = refineFlaggedTriangles(trinet, tri0, status)
    %
    % DESCRIPTION:
    %   Refine flagged triangles. Method depends on the flag degree, that
    %   is how many times the triangle was flagged, with
    %       - 1: Draw edge from split node to opposite corner
    %       - 2: Split the third edge as well, and refine by split in 4 as 
    %            in red step. Possibly flag neighbor triangle to new split.
    %       - 3: Split in 4 as in red step. 
    %
    % REQUIRED PARAMETERS:
    %   trinet  - TriNet (NetworkModel).
    %   tri0    - Triangulation before any refinement. 
    %   status  - Struct with info about refinements, flagging, edge splits
    %             etc. (see refineTriangles)
    %
    % RETURNS:
    %   trinet  - Updated TriNet
    %   status  - Updated status.
    
    % Iterate until no triangles are flagged
    while ~isempty(status.flagged)
        triID = status.flagged(1, 1);
        n = status.flagged(1, 2); % Number of refined neighbors (1-3)
        triangle = tri0(triID,:);
        triNodes = sort(triangle); % ascending order
        triEdges = [triNodes(1), triNodes(2);
                    triNodes(2), triNodes(3);
                    triNodes(1), triNodes(3)];
        split = status.splitNodes(all(status.splitNodes(:,1:2)==triEdges(1,:), 2) | ...
                                  all(status.splitNodes(:,1:2)==triEdges(2,:), 2) | ...
                                  all(status.splitNodes(:,1:2)==triEdges(3,:), 2), :);
        switch n
            case 1
                rededge = split(1:2);
                midnode = split(3);
                s = rededge(1);
                t = rededge(2);
                oppositeNode = triNodes(triNodes ~= s & ...
                                        triNodes ~= t);
                trinet = trinet.addEdge(oppositeNode, midnode);

                % Update tri
                newRow = find(ismember(trinet.tri, triangle, 'rows'));
                trinet.tri(newRow, :) = [];
                trinet.tri(end+1,:) = [s, oppositeNode, midnode];
                trinet.tri(end+1,:) = [t, oppositeNode, midnode];
                
            case 2
                % Flagged triangle has two "red" (previously splitted)
                % edges. Find node IDs from split:
                rededge1 = split(1,1:2); %  end nodes for first red edge
                midnode1 = split(1,3);   % split node for first red edge
                rededge2 = split(2,1:2); %  end nodes for second red edge
                midnode2 = split(2,3);   % split node for second red edge

                % The unsplit edge (end nodes):
                greenedge = triEdges(~all(triEdges == rededge1, 2) & ...
                                     ~all(triEdges == rededge2, 2),:); 

                % Need to determine which node of the green edge is on
                % which red edge. 
                % Impose that s is on rededge1 and t on rededge2.
                if ismember(greenedge(1), rededge1) && ismember(greenedge(2), rededge2)
                    s = greenedge(1);
                    t = greenedge(2);
                elseif ismember(greenedge(1), rededge2) && ismember(greenedge(2), rededge1)
                    s = greenedge(2);
                    t = greenedge(1);
                end

                % Draw edge from s to split node of second red edge, midnode2
                trinet = trinet.addEdge(s, midnode2);

                % Draw edge from midnode1 to midnode2
                trinet = trinet.addEdge(midnode1, midnode2);

                % Update tri
                oppositeNode = intersect(rededge1, rededge2);
                newRow = find(ismember(trinet.tri, triangle, 'rows'));
                trinet.tri(newRow, :) = [];
                trinet.tri(end+1,:) = [s, t, midnode2];
                trinet.tri(end+1,:) = [s, midnode1, midnode2];
                trinet.tri(end+1,:) = [midnode1, midnode2, oppositeNode];

            case 3
                midnode1 = split(all(split(:,1:2)==triEdges(1,:),2),3);
                midnode2 = split(all(split(:,1:2)==triEdges(2,:),2),3);
                midnode3 = split(all(split(:,1:2)==triEdges(3,:),2),3);
                trinet = trinet.addEdge(midnode1, midnode2);
                trinet = trinet.addEdge(midnode2, midnode3);
                trinet = trinet.addEdge(midnode1, midnode3);

                % Update tri
                newRow = find(ismember(trinet.tri, triangle, 'rows'));
                trinet.tri(newRow, :) = [];
                trinet.tri(end+1,:) = [triNodes(1), midnode1, midnode3];
                trinet.tri(end+1,:) = [triNodes(2), midnode1, midnode2];
                trinet.tri(end+1,:) = [triNodes(3), midnode2, midnode3];
                trinet.tri(end+1,:) = [midnode1, midnode2, midnode3];
        end
        status.flagged(status.flagged(:,1)==triID,:) = [];
    end
    status.flagged = [];
end

function status = flagTriangle(status, triID)
    if ~isempty(status.flagged) 
        if any(ismember(status.flagged(:,1), triID))
            % Flag opposite triangle for fixing later
            ix = find(status.flagged(:,1) == triID);
            status.flagged(ix, 2) = status.flagged(ix, 2) + 1; % Increase by one
        else
            status.flagged = [status.flagged; triID, 1];
        end
    else
        status.flagged = [triID, 1];
    end
end
