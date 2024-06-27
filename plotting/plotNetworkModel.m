function plotNetworkModel(netmod, varargin)
    % Plot network model.
    %
    % SYNOPSIS:
    %   plotNetworkModel(netmod, fig, 'pn1', pv1, ...)                      
    %   plotNetworkModel(netmod, 'pn1', pv1, ...)                      
    %
    % OPTIONAL PARAMETERS:
    %   See opt below.
    %

    if mod(nargin,2) == 0
        varargin = varargin(2:end);
        hold on
    else
        figure;
    end
    
    opt = struct('plotCoords',     true, ...  % if true, use physical coordinates
                 'nodesize',          4, ...
                 'linewidth',       0.5, ...  % line width for edges
                 'highlightwells', true, ...  % highlight wells with color and size
                 'wellcolor',       'r', ...  % well node color
                 'wellsize',          5, ...  % size of well nodes
                 'wellnamesize',     12, ...  % font size for well labels
                 'axis',          false, ...  % include axis
                 'label',       'wells', ...  % node labels 'wells', 'none' or 'numbers'
                 'voronoi',       false, ...  % include voronoi diagram
                 'flip',          false, ...  % switch x,y (rotate 90 degrees)
                 'pvScaled',        false);   % scale node sizes according to pore volume
    opt = merge_options(opt, varargin{:});

    network = netmod.network.G;
    grid = netmod.model.G;

    if opt.pvScaled
        pv = netmod.model.operators.pv;
        meanSize = opt.wellsize;
        fac = meanSize/mean(pv);  % mean pv will give normal size
        sizes = pv .* fac;
        opt.nodesize = sizes;
    end

    h = plot(network, 'MarkerSize', opt.nodesize, 'LineWidth', opt.linewidth);

    if opt.plotCoords && size(grid.cells.centroids,2) == 2
        if ~opt.flip
            h.XData = netmod.x;
            h.YData = netmod.y;
        else
            h.XData = netmod.y; % flips plot 90 degrees
            h.YData = -netmod.x;
        end
        view(2);
    elseif opt.plotCoords %3D case
        h.XData = netmod.x;
        h.YData = netmod.y;
        h.ZData = netmod.z;
    end
    axis equal

    % Highlight wells using different color and/or size.
    if opt.highlightwells
        if opt.pvScaled % do not change node size for wells
            highlight(h, netmod.wellNodes, ...
                    'NodeColor', opt.wellcolor, ...
                    'NodeFontSize', opt.wellnamesize)
        else
            highlight(h, netmod.wellNodes, ...
                    'NodeColor', opt.wellcolor, ...
                    'NodeFontSize', opt.wellnamesize, ...
                    'MarkerSize', opt.wellsize)
        end
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

    if opt.voronoi
        voronoiGrid = netmod.makeGrid();
        plotGrid(voronoiGrid, 'FaceColor', 'none', 'EdgeColor', [.7 .7 .7])
    end

    if ~opt.axis
        axis tight
        %axis off
        set(gca, 'xcolor', 'none')
        set(gca, 'ycolor', 'none')
    end
end
