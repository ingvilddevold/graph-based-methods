function plotNetworkModelCircular(netmod, varargin)
    % Plot network model using a circular layout.
    %
    % SYNOPSIS:
    %   plotNetworkModelCircular(netmod, fig, 'pn1', pv1, ...)                      
    %   plotNetworkModelCircular(netmod, 'pn1', pv1, ...)                      
    %
    % OPTIONAL PARAMETERS:
    %   See opt below.
    %

    if mod(nargin,2) == 0
        varargin = varargin(2:end);
    else
        figure;
    end

    opt = struct('nodesize',          4, ...
                 'highlightwells', true, ...
                 'wellcolor',       'r', ...
                 'wellsize',          5, ...
                 'wellnamesize',     12, ...
                 'label',        'wells');
    opt = merge_options(opt, varargin{:});

    network = netmod.network.G;
    
    h = plot(network, 'MarkerSize', opt.nodesize, ...
                      'layout', 'circle');

    % Highlight wells using different color and/or size.
    if opt.highlightwells
        highlight(h, netmod.wellNodes, ...
                     'NodeColor', opt.wellcolor, ...
                     'NodeFontSize', opt.wellnamesize, ...
                     'MarkerSize', opt.wellsize)
    end
    
    switch opt.label
        case 'wells'
            % Label well nodes only.
            for i = 1:netmod.numWells
                well = netmod.schedule.control(1).W(i);
                labelnode(h, well.cells, well.name);
            end
            % Set empty labels for non-well nodes.
            nonWellNodes = setdiff(1:network.numnodes, netmod.wellNodes);
            labelnode(h, nonWellNodes, repmat("", 1, numel(nonWellNodes)));
        case 'none'
            labelnode(h, 1:network.numnodes, repmat("", 1, network.numnodes))
        case 'numbers'
            labelnode(h, 1:network.numnodes, 1:network.numnodes);
        otherwise
            warning('Unsupported labeling option') 
    end
end
