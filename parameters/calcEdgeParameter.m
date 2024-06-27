
function v = calcEdgeParameter(network, edgeID, val, method)
    % Calculates parameter 'param' for node 'nodeID'
    % using specified method.
    switch method
        case 'avg_neighbors'
            v = edgeParamAvgNeighbor(network, edgeID, val);

        case 'avg_all'
            edges = 1:numEdges(netmod);
            edges = edges(edges ~= edgeID); % Exclude new edge
            v = mean(val(edges), 'omitnan');

        case 'sum_neighbors'
            v = edgeParamSumNeighbor(network, edgeID, val);
            
        case 'nan'
            v = nan;

        otherwise
            error('Unsupported fill-in method: %s', method)
    end
end    

%------------------------------------------------------------------
function v = edgeParamAvgNeighbor(network, edgeID, val)
    % Calculates parameter 'param' for edge 'edgeID' averaging
    % over the edges of the end nodes.
    [s,t] = network.getEndNodes(edgeID);
    sEdges = network.getEdgesFromNode(s); % Edges from s
    tEdges = network.getEdgesFromNode(t); % Edges from t
    edges = [sEdges; tEdges];         
    edges = unique(edges);
    edges = edges(edges~=edgeID);     % Remove edge of interest
    assert(~isempty(edges), 'End nodes have no other neighbors.')
    v = mean(val(edges), 'omitnan');
    % Prevent break-down when new edge only has new edges as
    % neighbors (giving val = NaN) by using global average:
    if isnan(v)
        edges = 1:network.numEdges;
        edges = edges(edges ~= edgeID); % Exclude new edge
        v = mean(val(edges), 'omitnan');
    end
end

%------------------------------------------------------------------
function val = edgeParamSumNeighbor(network, edgeID, val)
    % Calculates parameter 'param' for edge 'edgeID' summing
    % over the edges of the end nodes.
    [s,t] = network.getEndNodes(edgeID);
    sEdges = network.getEdgesFromNode(s);
    tEdges = network.getEdgesFromNode(t);
    edges = [sEdges; tEdges];
    edges = unique(edges);
    edges = edges(edges~=edgeID);
    assert(~isempty(edges), 'End nodes have no other neighbors.')
    val = sum(val(edges), 'omitnan');
end
