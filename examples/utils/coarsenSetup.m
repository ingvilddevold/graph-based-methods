function cSetup = coarsenSetup(setup, celldim)
    % Make a coarse simulation model from a fine.
    %
    % SYNOPSIS:
    %   cSetup  = coarseModel(setup, [6 6 1])    6x6x1 Cartesian grid
    %
    % REQUIRED PARAMETERS:
    %   setup      - A SimulatorSetup with a model, schedule and state0.
    %   celldim    - Vector, length 3, giving the number of 
    %                cells in each direction.
    %  
    % RETURNS:
    %   cSetup   -  New SimulatorSetup for the coarsened model.
    %

    % Fine-scale model, schedule, state0.
    fModel     = setup.model;
    fSchedule  = setup.schedule;
    fState0    = setup.state0;
    
    % Make coarse grid.
    G = fModel.G;
    nmin = min(G.nodes.coords)-1e-5;
    nmax = max(G.nodes.coords)+1e-5;
    cG = cartGrid(celldim, nmax-nmin);
    cG.nodes.coords = cG.nodes.coords + nmin;
    cG = computeGeometry(cG);
    ix = findEnclosingCell(cG, G.cells.centroids);
    cno = accumarray(ix,1,[cG.cells.num 1]);
    cG = extractSubgrid(cG, cno>0);
    
    crock = makeRock(cG, mean(fModel.rock.perm), ...
        mean(fModel.rock.poro));
    cModel =  fModel;
    cModel.G = cG;
    cModel.rock = crock;
    cModel.FlowDiscretization = [];
    cModel.AutoDiffBackend = AutoDiffBackend();
    cModel = cModel.setupOperators(cG, cModel.rock);
    cModel.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

    % rescale pore volumes
    cModel.operators.pv = cModel.operators.pv * sum(setup.model.operators.pv)/sum(cModel.operators.pv);
    
    % Set the initial state and the simulation schedule
    cState0   = initResSol(cG, mean(fState0.pressure), mean(fState0.s));
    cSchedule = fSchedule;
    for n=1:numel(cSchedule.control)
        cW = [];
        for i=1:numel(cSchedule.control(n).W)
            oW = cSchedule.control(n).W(i);
            pt = fModel.G.cells.centroids(oW.cells,:);
            cellIx = unique(findEnclosingCell(cG, pt));
            cW = addWell(cW, cG, crock, cellIx(cellIx>0), 'Name', oW.name, ...
                'Type', oW.type, 'Val', oW.val, 'Radius', mean(oW.r), ...
                'dir','z', 'compi', oW.compi, 'Sign', oW.sign,...
                'status', oW.status, 'lims', oW.lims);
        end
        cSchedule.control(n).W = cW;
    end
    cSetup = struct('model', cModel, 'schedule', cSchedule, 'state0', cState0);
end

