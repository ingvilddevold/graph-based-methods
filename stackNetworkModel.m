function netmod = stackNetworkModel(netmod2d, heights, params)
    % Create a 2.5D NetworkModel by stacking 2D models
    %
    % SYNOPSIS:
    %   netmod = stackNetworkModel(netmod2d, heights)
    %
    % DESCRIPTION: 
    %   The 2D network is repeated in a stack with layers at specified
    %   heights, and vertical connections between the layers.
    %   Parameters are filled in, the well structure updated to match
    %   new number of perforations and state0 expanded to match new 
    %   number of nodes/cells.
    %
    % REQUIRES PARAMETERS:
    %   netmod2d - Two-dimensional NetworkModel
    %   heights  - Array of heights for the layers
    %
    % RETURNS:
    %   netmod  - 2.5D NetworkModel. Ready to simulate
    %

    netmod = netmod2d;
    nLayers = numel(heights);

    [network, newnodes, newedges] = stackNetwork(netmod2d.network, heights);

    netmod.network = network;

    %%% UPDATE MODEL
    netmod.model = netmod2d.model;
  
    netmod.model.G.cells.num = network.numNodes;
    netmod.model.G.cells.centroids = netmod.getNodeCoordinates();

    netmod.model = netmod.model.removeStateFunctionGroupings();

    netmod = fillinParameters(netmod, params, nLayers);

    netmod.model = netmod.model.setupOperators(netmod.model.G, ...
        netmod.model.rock, ...
        'trans', netmod.model.operators.T, ...
        'neighbors', network.G.Edges.EndNodes, ...
        'porv', netmod.model.operators.pv);

    %%% UPDATE SCHEDULE
    % Create new well structure with correct number of perforations
    schedule = netmod.schedule;
    for n = 1:numel(schedule.control)
        W = [];
        for i = 1:numel(schedule.control(n).W)
            oW = schedule.control(n).W(i);
            if numel(oW.r) < 2
                % use same radius at all perforations
                r = repmat(oW.r, nLayers, 1);
            else
                r = oW.r(1:nLayers);
            end
            wellCells = [oW.cells];
            for l = 2:nLayers
                wellCells = [wellCells; oW.cells + netmod2d.network.numNodes*(l-1)];
            end
            W = addWell(W, netmod.model.G, netmod.model.rock, wellCells, ...
                'Type',        oW.type, ...
                'Val',         oW.val, ...
                'Radius',      r, ...
                'Name',        oW.name, ... 
                'sign',        oW.sign, ...
                'lims',        oW.lims, ...
                'status',      oW.status, ...
                'refDepth',    oW.refDepth, ...
                'calcReprRad', false, ...
                'compi',       oW.compi, ...
                'WI',          repmat(oW.WI, nLayers, 1));
        end
        schedule.control(n).W = W;
    end
    netmod.schedule = schedule;

    %%% UPDATE STATE0
    netmod.state0 = netmod.expandState0();
end

function [network, newnodes, newedges] = stackNetwork(network, heights)
    % Go from 2D to 2.5D network.
    % 
    % SYNOPSIS: model3D = stack2Dmodel(model2D, heights)
    %
    % DESCRIPTION:
    %   Repeats the two-dimensional network at each
    %   z level / height in heights. 
    %   Adds vertical edges between the layers.
    %
    % REQUIRED PARAMETERS:
    %   network - Two-dimensional BaseGraph.
    %   heights - Array of heights / z coordinates for each layer.
    %
    % RETURNS: 
    %   network  - 2.5D BaseGraph.
    %   newnodes - nodeIDs for added nodes. (for filling in parameters)
    %   newedges - edgeIDs for added edges. 

    newnodes = [];
    newedges = [];

    N = network.numNodes;
    E = network.numEdges;
    
    % let the 2d model have the first height
    network.G.Nodes.z = repmat(heights(1), N, 1); 

    for i = 2:numel(heights) 
        zi = heights(i);
        
        % Add nodes and vertical connections
        for n = 1:N
            [network, nodeID]= network.addNode();

            newnodes = [newnodes; nodeID];

            % Set coordinates:
            network.G.Nodes.x(nodeID) = network.G.Nodes.x(n);
            network.G.Nodes.y(nodeID) = network.G.Nodes.y(n);
            network.G.Nodes.z(nodeID) = zi;
        end

        % Add horizontal connections
        for e = 1:E
            % For each edge in the original layer, add
            % the corresponding edge in the new layer.
            edgeToCopy = network.G.Edges.EndNodes(e,:);
            s = edgeToCopy(1) + N*(i-1);
            t = edgeToCopy(2) + N*(i-1);
            [network, edgeID] = network.addEdge(s,t);

            newedges = [newedges; edgeID];
        end
    end

    for i = 2:numel(heights)    
        % Add vertical connections
        for n = 1:N
            [network, edgeID] = network.addEdge(n, n+N);
            newedges = [newedges; edgeID];
        end
    end  
end

function netmod = fillinParameters(netmod, params, nLayers)
    % Fill in missing parameter values in new layers
    for i = 1:numel(params)
        val = netmod.getParameterValue(params{i});
        
        switch params{i}.name
            case 'porevolume'
                % redistribute vertically
                val = repmat(val/nLayers, nLayers, 1);
            case 'transmissibility'
                % use the same horizontal transmissibilities
                % use a mean for the new vertical connections
                % (vertical connections are last in edges)
                val = repmat(val, nLayers, 1);
                numNodesOrig = netmod.network.numNodes/nLayers;
                numVertEdges = numNodesOrig * (nLayers-1);
                meanT = mean(val);
                val = [val; repmat(meanT, numVertEdges, 1)];
        end
        netmod = netmod.setParameterValue(params{i}, val);
    end
    
    if isfield(netmod.model.rock, 'krscale')
        
        netmod.model.rock.krscale.drainage.ow = repmat( ...
                        netmod.model.rock.krscale.drainage.ow, nLayers, 1);
        netmod.model.rock.krscale.drainage.w = repmat( ...
                        netmod.model.rock.krscale.drainage.w, nLayers, 1);
        netmod.model.rock.krscale.drainage.og = repmat( ...
                        netmod.model.rock.krscale.drainage.og, nLayers, 1);
        netmod.model.rock.krscale.drainage.g = repmat( ...
                        netmod.model.rock.krscale.drainage.g, nLayers, 1);
    end
 end
