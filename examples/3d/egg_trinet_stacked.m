mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr

%% Set up Egg model
egg = TestCase('egg_wo');

problem = egg.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);
%egg.plot(statesFine)
%plotWellSols(wellSolsFine);

%% Need a 2D model to construct graph. Each well should only have one cell
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [25 25 1]);

%%
pts = cProblem.SimulatorSetup.model.fluid.krPts;
cProblem.SimulatorSetup.model = imposeRelpermScaling(cProblem.SimulatorSetup.model, ...
    'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
    'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
    'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
cProblem.SimulatorSetup.model = cProblem.SimulatorSetup.model.setupOperators( ...
    cProblem.SimulatorSetup.model.G, cProblem.SimulatorSetup.model.rock);

%% 2D TriNet
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0, ...
                'thres', 10);

%% Set up parameters and objective function

config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'porevolume',       1,   'linear',    [],       [],     [],        [0.01 10],   'node'
    'conntrans',        1,   'log',       [],       [],     [],        [0.01 100],  'well'
    'transmissibility', 1,   'log'        [],       [],     [],        [0.01 100],  'edge'
    'swl',              1,   'linear',    [0 .2],       [],     [],        [],      'node'
    'swcr',             1,   'linear',    [0 .3],       [],     [],        [],      'node'
    'swu',              1,   'linear',    [.7 1],       [],     [],        [],      'node'
    'sowcr',            1,   'linear',    [0 .3],       [],     [],        [],      'node'
    'krw',              1,   'linear',    [.5 2],       [],     [],        [],      'node'
    'kro',              1,   'linear',    [.5 2],       [],     [],        [],      'node'}; 

params = setupParameters(trinet, config);

%% Define objective function
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);


%% Stacked TriNet

trinet.model.operators.gdz = trinet.model.getGravityGradient();

top = min(egg.model.G.nodes.coords(:,3));
bottom = max(egg.model.G.nodes.coords(:,3));
trinet_stack = stackNetworkModel(trinet, [top, bottom], params);

params = setupParameters(trinet_stack, config);

trinet_stack.model.operators.gdz = trinet_stack.model.getGravityGradient();

%% gdz

% Set box lims based on top and bottom z coordinate
top = min(egg.model.G.cells.centroids(:,3));
bottom = max(egg.model.G.cells.centroids(:,3));
ulim = (bottom-top)*15;

config_gdz = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'gdz',              1,   'linear',  [-ulim, ulim],   [],     [],    [], 'edge'}; 
params = addNetworkModelParameter(params, trinet_stack, ...
        'name',    config_gdz{1,1}, 'scaling', config_gdz{1,3}, ...
        'boxLims', config_gdz{1,4}, 'lumping', config_gdz{1,5}, ...
        'subset',  config_gdz{1,6}, 'relativeLimits',config_gdz{1,7}, ...
        'uniformLimits', false, ...
        'belongsTo', 'model', ...
        'location', {'operators' 'gdz'}, ...
        'mapTo', 'edge');

%%
trinet_stack.model = trinet_stack.reimposeRelpermScaling();

%%
trinet_stack.params = params;
problem = trinet_stack.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name', 'egg_trinet_stacked', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_stack_tuned, p, h] = optimizeNetworkModel(trinet_stack, p0, 'maxIt', 40);
