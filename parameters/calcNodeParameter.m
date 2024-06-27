function v = calcNodeParameter(network, nodeID, val, method)
    % Calculate node parameter value
    %
    % REQUIRED PARAMETERS:
    %   network - BaseGraph.
    %   nodeID  - Node missing parameter value
    %   val     - Array of parameter values
    %   method  - Fill-in strategy

    switch method
        case 'avg_neighbors'
            nodes = network.getNeighbors(nodeID);
            if isempty(nodes)
                warning('Node has no neighbors. Using avg_all instead.')
                nodes = 1:numNodes(network);
                nodes = nodes(nodes~=nodeID); 
            end
            v = mean(val(nodes), 'omitnan');

            if isnan(v)
                warning(['Neighbor averaging gave nan parameter value. Using' ...
                    'global average instead.'])
                v = calcNodeParameter(network, nodeID, val, 'avg_all');
            end
    
        case 'avg_all'
            nodes = 1:numNodes(network);
            nodes = nodes(nodes~=nodeID); 
            v = mean(val(nodes), 'omitnan');
        
        case 'sum_neighbors'
            neighIDs = obj.network.getNeighbors(nodeID);
            v = sum(val(neighIDs), 'omitnan');
        
        case 'nan'
            v = nan;
    
        otherwise
            error('Unsupported fill-in method: %s', method)
    end
end


