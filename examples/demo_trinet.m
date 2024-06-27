mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr

%% Set up demo water-oil model

test = TestCase('diagnostics_2d_wo', 'lognormal', true, 'barriers', false, 'pvi', 1.5);
%test.plot('plotWells', true)

problem = test.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);
%test.plot(statesFine)
%plotWellSols(wellSolsFine);

%% Set up TriNet model
trinet = TriNet(test.model, test.schedule, test.state0, 'useDistmesh', true);

%% Set up parameters and objective function

problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num', 1);

config = {
    ...%name      include   scaling  boxlims lumping subset relativeLimits  mapTo
    'porevolume',       1, 'linear',   [],     [],    [],     [0.01 10],   'node' 
    'conntrans',        1, 'log',      [],     [],    [],     [0.01 100],  'well'
    'transmissibility', 1, 'log'       [],     [],    [],     [0.01 100],  'edge'}; 
params = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end     % include = 0
    params = addNetworkModelParameter(params, problem.SimulatorSetup, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'lumping', config{k,5}, ...
        'subset',  config{k,6}, 'relativeLimits',config{k,7}, ...
        'uniformLimits', false);
end

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

trinet.params = params;

%% Triangle refinement

p0 = OptimizationProblem(samples, ...
            'parameters',        params, ...
            'name', 'diagnostics_wo_trinet_ref', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet, p, h] = optimizeNetworkModelTopology(trinet, p0, config, ...
                'refineFrom', 'triangles', 'selectWith', 'fraction', ...
                'maxFrac', 0.10, 'sensMap', 'mean', ...
                'maxIt', 200, 'outerMaxIt', 20, 'totalMaxIt', 1000, ...
                'outerObjTol', 1e-4, 'resTolAbs', 1e-9, ...
                'name', "TriNet (dynamic)", 'resChangeTolRel', 0.01);

%% Static tuning

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',  'diagnostics_wo_trinet_static', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet, p, h] = optimizeNetworkModel(trinet, p0, 'maxIt', 40);
