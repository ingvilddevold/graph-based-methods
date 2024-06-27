function [tri, G, wellNodes, state0] = delaunayGraph(model, schedule, state0, varargin)
    % Make a simulation graph using Delaunay triangulation.
    %
    % SYNOPSIS:
    %   [tri, G, wellNodes, state0] = delaunayGraph(model, schedule, state0)
    %   
    % DESCRIPTION:
    %   Given a regular simulation problem with a model, schedule and 
    %   state0, construct a graph by a Delaunay triangulation of 
    %   the wells and the outer boundary.
    %
    % REQUIRED PARAMETERS:
    %   model    - Fluid model.
    %   schedule - Simulation schedule incl. wells and time steps.
    %   state0   - Initial state.
    %
    % OPTIONAL PARAMETERS:
    %   See delaunayTriangulate.
    %
    % RETURNS:
    %   tri         - Triangulation.  Each row contains three point indices.
    %   G           - Graph based on triangulation. 
    %   wellNodes   - Array of nodes corresponding to well cells.
    %   state0      - Initial state matching graph dimensions. 
    
    [tri, coords, wellNodes] = delaunayTriangulate(model, schedule, varargin{:});
    G = triToGraph(tri);
    dim = size(model.G.cells.centroids, 2);
    coordnames = ['x', 'y', 'z'];
    for i = 1:dim
        G.Nodes.(coordnames(i)) = coords(:,i);
    end
    state0 = makeState0(G, state0);
end

%--------------------------------------------------------------------------
function [tri, coords, wellNodes] = delaunayTriangulate(model, schedule, varargin)
    % Delaunay triangulation of wells and outer boundary.
    % 
    % SYNOPSIS:
    %   [tri, coords, wellNodes] = delaunayTriangulate(model, schedule)
    %
    % DESCRIPTION:
    %   Constructs a Delaunay triangulation using the boundary points /
    %   convex hull of the domain and the wells. This triangulation can 
    %   then be used to construct a graph.
    %   If 'useDistmesh' is true, DistMesh is used before the
    %   Delaunay triangulation. 
    %   When defining the outer boundary points using the
    %   convex hull, points closer than ´thres´ will be merged.
    %
    % REQUIRED PARAMETERS:
    %   model    - Fluid model.
    %   schedule - Simulation schedule incl. wells and time steps.
    %
    % OPTIONAL PARAMETERS:
    %   'thres'       - Threshold for distance between points in convex hull.
    %   'useDistmesh' - Boolean. Whether to adjust points with DistMesh.
    %   'edgeFac'     - Max edge length for DistMesh
    %   'useNodes'    - Boolean. Use nodes instead of centroids.
    %   'fixBoundary' – Boolean. Fix both wells and boundary points in DistMesh.
    %   'hMin'        - DistMesh minimum relative edge length
    %   'hMax'        - DistMesh maximum relative edge length
    %
    % RETURNS:  
    %   tri         - Triangulation.  Each row contains three point indices.
    %   coords      - Coordinates of the triangulation points.
    %   wellNodes   - Indices for the nodes corresponding to cells
    %                 containing wells.

    opt = struct('thres', 100, ...
                 'useDistmesh', true, ...
                 'adaptToWells', false, ...
                 'edgeFac', 0.25, ...
                 'useNodes', false, ...
                 'fixBoundary', false, ...
                 'hMin', 0.075, ...
                 'hMax', 0.5);

    opt = merge_options(opt, varargin{:});
   
    G = model.G;
    dim = size(G.cells.centroids, 2);

    %----------------------------------------------------------------------
    % Find outer boundary using convex hull.
    % Using cell corner nodes if the grid has them, otherwise centroids.
    if opt.useNodes && isfield(G, 'nodes')
        xb = G.nodes.coords;
        % Merge the points with equal x,y coordinates to a point
        % in the middle, s.t. z = (z1 + z2)/2
        % Will then get z coordinates consistent with those for wells.
        xb_sorted = sortrows(xb, [1 2]);
        if dim==3
            xb_lower = xb_sorted(1:2:end, :);
            xb_upper = xb_sorted(2:2:end, :);
            xb = [xb_lower(:,1), xb_lower(:,2), (xb_upper(:,3)+xb_lower(:,3))/2]; 
        elseif dim==2 %2D case
            xb = xb_sorted;
        end
    else
        % Use centroids instead of nodes. 
        xb = G.cells.centroids;
    end
    
    ch = convhull(xb(:,1:2), 'simplify', true); % Indices of nodes in 2D convex hull
    xb = xb(ch,:);                              % Coordinates for nodes in convex hull 
    % Note: first and last in xb are the same node

    %----------------------------------------------------------------------
    % Merge the boundary points that are closer than opt.thres, by
    % adding one point in middle instead of the two.
    % convhull gives the points in counterclockwise order,
    % so we only check subsequent points.
    dist = sqrt( sum( abs( diff( xb(:,1:2)) ).^2, 2 ) ); % distance in xy-plane
    merge = dist < opt.thres;
    
    xb2 = zeros(size(xb,1) - 1 - sum(merge), dim);
    j = 1;
    for i = 1:numel(merge)
        if i > 1 && merge(i-1) % point was merged in previous iteration
            continue;
        elseif merge(i)
            if i == numel(merge) % last and first are too close
                xb2(1,:) = (xb(i,:) + xb(1,:))/2;
            else
                xb2(j,:) = (xb(i,:) + xb(i+1,:))/2;
                j = j+1;
            end
        else    % accept point
            xb2(j,:) = xb(i,:); 
            j = j+1;
        end
    end
    xb = xb2;

    %----------------------------------------------------------------------
    % Define the points to triangulate around
    
    % Wells
    wc = arrayfun(@(W) W.cells(1), schedule.control(1).W);  % cell indices
    xw = unique(G.cells.centroids(wc,:), 'rows', 'stable'); % coordinates
    
    % Concatenate x and y coordinates
    x = [xw(:,1:2); xb(:,1:2)];
    
    %x = xw(:,1:2); % SPE10 hack

    %----------------------------------------------------------------------
    % Triangulate
    if opt.useDistmesh
        bbox = [min(xb); max(xb)];              % Bounding box
        ldiag = sqrt(sum(diff(bbox,1,1).^2,2)); % Bounding box diagonal length
        h0 = ldiag * opt.edgeFac;               % Initial edge length
        
        if opt.adaptToWells
            % Distance from x to nearest well
            d  = @(x) minPdist2(x,xw(:,1:2));
            % Function that goes to 0 near xw and 1 "sufficiently far away". Can be
            % changed to any suitable function
            t  = @(x) min(exp(d(x)/(0.25*ldiag)) - 1, 1);
            % Relative edge length
            hMin = opt.hMin;                         % Minimum relative edge length
            hMax = opt.hMax;                         % Maximum relative edge length
            fh = @(x, varargin) t(x).*hMax + (1-t(x)).*hMin;
        else
            fh = @(x, varargin) ones(size(x,1),1);
        end

        if opt.fixBoundary
            xFix = x;
        else
            xFix = xw(:,1:2);
        end

        [x, ~, sorting] = distmesh2d(...
            @dpoly, ...                         % Distance function
            fh, ...                             % Scaled edge length (1 = uniform)
            h0, ...                             % Initial edge length
            bbox, ...                           % Boundary
            xFix, ...                           % Fixed points
            false, ...                          %
            [xb; xb(1,:)]);                     % Boundary
        
        % DistMesh changes the order of all sites - undo this sorting.
        isWell  = false(max(sorting),1); isWell(1:size(xw,1))= true;
        isWell  = isWell(sorting);
        [~,ixw] = sort(sorting(isWell));
        xw = x(isWell,:); xw = xw(ixw,:);
        xr = x(~isWell,:);
        % Triangulation points
        x = [xw; xr];

        coords = x;
        % Workaround for equal z coordinates:
        if dim==3
            coords(:,3) = repmat(mean(xb(:,3)), size(coords,1),1);
        end
    else
        coords = zeros(size(xw,1)+size(xb,1), dim);
        for i = 1:dim
            coords(:,i) = [xw(:,i); xb(:,i)];
        end
    end
    % Delaunay triangulation
    tri = delaunay(x);

    % if delaunay uses points in same order as in x, 
    % should work with
    wellNodes = 1:size(xw,1);
end

%--------------------------------------------------------------------------
function G = triToGraph(tri)
    % Make a MATLAB graph from a triangulation.
    %
    % SYNOPSIS:
    %   G = triToGraph(tri)
    %
    % DESCRIPTION:
    %   Converts a tri matrix to an edge matrix and makes a graph.
    %   Each row in tri represents a triangle with the three columns 
    %   being its three corner points.
    %   Each point is mapped to a node in the graph, and the sides of 
    %   the triangles are mapped to edges. 
    %
    % REQUIRED PARAMETERS:
    %   tri - Triangulation. Array of size #Triangles x 3
    %
    % RETURNS:
    %   G - graph with the same geometry as the triangulation. 
   
    % Extract edge info
    N1 = tri(:,1:2);
    N2 = tri(:,2:3);
    N3 = tri(:,[1,3]);
    N = [N1; N2; N3];

    % Impose convention of smallest node first
    % to avoid duplicated edges (opposite directions)
    % (There is probably a much better way to do this)
    for r = 1:size(N,1)
        if N(r,1) > N(r,2)
            s = N(r,2);
            N(r,2) = N(r,1);
            N(r,1) = s;
        end
    end
    % Remove duplicates
    N = unique(N, 'rows');

    G = graph(N(:,1), N(:,2));
end

%--------------------------------------------------------------------------
function tristate0 = makeState0(G, state0)
    % Make a state0 matching the dimensions of the graph G.
    %
    % SYNOPSIS:
    %   tristate0 = makeState0(G, state0)
    %
    % DESCRIPTION:
    %   Make a new initial state setting all values to the means 
    %   of the original initial state.
    %
    % REQUIRED PARAMETERS:
    %   G           - Graph (constructed from triangulation).
    %   state0      - The initial state of the original simulation problem. 
    %
    % RETURNS:
    %   triState0   - New initial state matching the graph G. 

    nn = size(G.Nodes, 1); % number of nodes
    ne = size(G.Edges, 1); % number of edges
    stateDim = struct('pressure', nn, ...
                      'flux',     ne, ...
                      's',        nn);
    tristate0 = struct;
    fn = fieldnames(state0);
    for i = 1:numel(fn)
        tristate0.(fn{i}) = repmat(mean(state0.(fn{i})), [stateDim.(fn{i}), 1]);
    end
end
