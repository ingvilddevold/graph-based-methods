function [netmod, p, histAll] = optimizeNetworkModelTopology(netmod, p, paramConfig, varargin)    
    % Main graph optimization routine.
    %
    % SYNOPSIS:
    %   [p, h] = optimizeNetworkModelTopology(netmod, p, paramConfig, 'pn1', pv1, ...)
    %
    % DESCRIPTION:
    %   Attempts to optimize both the parameter values and the layout of
    %   the graph itself. (Graph refinement is only supported for TriNet)
    %   Parameters are tuned with the Levenberg-Marquardt 
    %   algorithm (unitBoxLM.m) for each graph.
    %   After parameter tuning, triangles are selected for refinement based
    %   on the standard deviation of the sensitivities of the mismatch 
    %   function wrt. the parameters. For details on the refinement itself,
    %   see refineGraph.
    %
    %   We refer to the graph refinement as the outer loop and the
    %   parameter tuning as the inner loop.
    %   
    % REQUIRED PARAMETERS:
    %   netmod      - NetworkModel
    %   p           - Optimization problem.
    %   paramConfig - Parameter config cell array.
    %
    % RETURNS:
    %   netmod      - NetworkModel
    %   p           - Optimization problem (last).
    %   histAll     - Struct containing mismatch and projected gradient
    %                 values in all iterations.
    %
    % OPTIONAL PARAMETERS:
    %   See opt, refineOpt below. 
    %   In addition any arguments supported by optimize and unitBoxLM, e.g.
    %   'plotEvolution' - Set to false to avoid progress plots in each
    %                     parameter tuning.
    %

    opt = struct('outerMaxIt',       5, ...  % max iterations in outer loop
                 'totalMaxIt',     200, ...  % max total inner iterations    
                 'outerGradTol',  1e-6, ...  % gradient tolerance in outer loop
                 'outerObjTol',   1e-6, ...  % objective tolerance in outer loop
                 'plotGraphs',    true, ...  % plot graph in each outer iteration
                 'saveGraphs',    true, ...  % save graph in each outer iteration
                 'saveWellSols',  true, ...  % save well sols after each outer iteration
                 'plotWellSols', false, ...  % plot well curves after each outer iteration
                 'wellSolsFine',    [], ...  % fine-scale well sol for plot^        
                 'plotHist',      true, ...  % plot mismatch/gradient history across all outer iterations
                 'name',            '', ...  % name of refinement case for plot legend
                 'refineGraph', true);         
    
    % Options for refinement
    refineOpt = struct('refineFrom', 'triangles', ...    
                       'sensMap',         'mean', ...
                       'selectWith',  'fraction', ...
                       'maxSens',            0.5, ...
                       'maxFrac',            0.1);

    [opt, varargin] = merge_options(opt, varargin{:});
    [refineOpt, optimConfig] = merge_options(refineOpt, varargin{:});
    
    it = 0;     % outer iterations / graph refinements
    itAll = 0;  % total iterations                  
    n = opt.outerMaxIt+1;
    h = struct('val', nan(1,n), 'pg', nan(1,n));
    histAll = struct('val', nan(1,n), 'pg', nan(1,n), ...
                     'its', nan(1,n), 'name', opt.name);
    dir = p.mainDirectory;
    if opt.plotHist
        fig = figure;
    end

    while ~converged(it, itAll, h, opt)
        it = it + 1;
        
        % Update graph
        %------------------------------------------------------------------
        if it > 1 && opt.refineGraph
            [netmod, rescaleNodes] = netmod.refineGraph(refineOpt);
            netmod = netmod.fillinNanParameters(rescaleNodes);
            netmod.state0 = netmod.expandState0();
            netmod.model = netmod.updateModel();
        end

        % Plot graph
        if opt.plotGraphs
            netmod.plot();
        end
        %------------------------------------------------------------------

        % Tune graph
        %------------------------------------------------------------------
        % Set up new optimization problem p
        params = makeParams(netmod, paramConfig);
        netmod.params = params;

        problem = netmod.getPackedSimulationProblem();
        samples = struct('problem', {{problem}}, ...
                         'num',     1);
        
        hash = obj2hash(problem);
        directory = fullfile(dir, hash);
        p = OptimizationProblem(samples, ...
                                'parameters', params, ...
                                'name',       p.name, ...
                                'directory',  directory, ...
                                'objective',  p.objective, ...
                                'setupType',  'simulation', ...
                                'verboseSimulation', false, ...
                                'solverFunOptions',  {'scalarObjective', false});
        % Run optimization
        [netmod, p, hist] = optimizeNetworkModel(netmod, p, optimConfig{:});
        p = calcParamSensitivities(p);

        % Update graph parameter and sensitivitiy values
        netmod.params = p.parameters;
        
        % Write graph to disk
        if opt.saveGraphs
            dest = fullfile(dir, 'graphs');
            mkdir(dest)
            filename = sprintf('graph_%d.mat', it);
            save(fullfile(dest, filename), 'netmod')
        end

        [h.val(it), h.pg(it)] = deal(hist.val(end), hist.pg(end));
        nval = numel(hist.val);
        [histAll.val(itAll+1:itAll+nval), histAll.pg(itAll+1:itAll+nval)] = deal(hist.val, hist.pg);
        histAll.its(it) = nval;
        itAll = itAll + nval;

        if opt.plotHist
            fig = plotHist(fig, histAll, 'legend', false);
        end

        if opt.saveWellSols
            [wellSols, ~] = getPackedSimulatorOutput(p.samples.problem{1});
            dest = fullfile(dir, 'wellsols');
            mkdir(dest)
            filename = sprintf('wellsols_%d.mat', it);
            save(fullfile(dest, filename), 'wellSols')
        end

        if opt.plotWellSols
            simulatePackedProblem(p.samples.problem{1});
            [wellSolsNew, ~] = getPackedSimulatorOutput(p.samples.problem{1});
            plotWellSols({opt.wellSolsFine, wellSolsNew}, 'datasetnames', {'true', 'predicted'});

            dest = fullfile(dir, 'wellcurves');
            mkdir(dest)
            filename = sprintf('wellcurve_%d.fig', it);
            savefig(fullfile(dest, filename))
        end
        % Save history to output directory
        save(fullfile(dir, 'history.mat'), 'histAll')
    end
end

%--------------------------------------------------------------------------
function flag = converged(it, itAll, h, opt)
    flag = false;
    if it > 0
        flags = [it >= opt.outerMaxIt, ...
                 itAll >= opt.totalMaxIt, ...
                 h.pg(it) < opt.outerGradTol, ...
                 h.val(it) < opt.outerObjTol];
        flag = any(flags);
    end 
    if flag
        ll = 65;
        prnt = @(s)fprintf('| %s%s |\n', s, repmat(' ', [1, ll-length(s)-4]));
        fprintf('%s\n', repmat('-', [1, ll]));
        prnt('Graph optimization finished:');
        switch find(flags,1)
            case 1
                s = sprintf('Reached maximum number of outer iterations/ graph refinements (%d)', it);
            case 2
                s = sprintf('Reached maximum number of total iterations (%d)', itAll);
            case 3
                s = sprintf('Norm of projected gradient below tolerance (%7.2e < %7.2e)', h.pg(it), opt.outerGradTol);
            case 4
                s = sprintf('Mismatch below tolerance (%7.2e < %7.2e)', h.val(it), opt.outerObjTol);
        end
        prnt(s)
        fprintf('%s\n', repmat('-', [1, ll]));
    end
end

%--------------------------------------------------------------------------
function params = makeParams(setup, config)
    params = [];
    for k = 1:size(config,1)
        if config{k, 2} == 0, continue, end     % include = 0
        params = addNetworkModelParameter(params, setup, ...
            'name',    config{k,1}, 'scaling', config{k,3}, ...
            'boxLims', config{k,4}, 'lumping', config{k,5}, ...
            'subset',  config{k,6}, 'relativeLimits',config{k,7}, ...
            'uniformLimits', false, 'mapTo', config{k, 8});
    end
end

%--------------------------------------------------------------------------
function p = calcParamSensitivities(p)
    % Based on updateModelParameters
    % Maps sensitivities from matrix sens to 
    % each parameter in a cell array.
    pinit = getScaledParameterVector(p.samples.problem{1}.SimulatorSetup, p.parameters);
    [r, J] = p.getScaledObjective(pinit);
    sensMatrix = std(J)';
    nparam = cellfun(@(x)x.nParam, p.parameters);
    u = mat2cell(sensMatrix, nparam, 1);
    for k = 1:numel(p.parameters)
        p.parameters{k}.sens = u{k};
    end
end
