mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr

%% Set up SAIGUP case

saigup = TestCase('saigup_wo');

% uniform saturation
saigup.state0.s = saigup.state0.s.*0 + [0,1];

problem = saigup.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);

%% Need a 2D model to construct graph. Each well should only have one cell
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [25 25 1]);

%% 2D TriNet
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0);
 
%% Set up parameters and objective function

config = {
    ...%name      include   scaling  boxlims lumping subset relativeLimits  mapTo
    'porevolume',       1, 'linear',   [],     [],    [],     [0.01 10],   'node' 
    'conntrans',        1, 'log',      [],     [],    [],     [0.01 100],  'well'
    'transmissibility', 1, 'log'       [],     [],    [],     [0.01 100],  'edge'}; 

params = setupParameters(trinet, config);


%% Define objective function
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% 2D TriNet

problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',    'saigup_trinet', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_tuned, p, h] = optimizeNetworkModel(trinet, p0, 'maxIt', 40);

%% Finer TriNet
trinet_fine = TriNet(cProblem.SimulatorSetup.model, ...
                     cProblem.SimulatorSetup.schedule, ...
                     cProblem.SimulatorSetup.state0, ...
                     'useDistmesh', true, 'edgeFac', 0.1);

params = setupParameters(trinet_fine, config);
trinet_fine.params = params;

problem = trinet_fine.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name', 'saigup_trinet_edge0_1', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_tuned, p, h] = optimizeNetworkModel(trinet_fine, p0, 'maxIt', 40);

%% Stacked TriNet

trinet.model.operators.gdz = trinet.model.getGravityGradient();
trinet_stack = stackNetworkModel(trinet, [2200, 2500], params);
trinet_stack.model.operators.gdz = trinet_stack.model.getGravityGradient();

params = setupParameters(trinet_stack, config);

%% gdz

% Set box lims based on top and bottom z coordinate
top = min(saigup.model.G.cells.centroids(:,3));
bottom = max(saigup.model.G.cells.centroids(:,3));
ulim = (bottom-top)*10;

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
trinet_stack.params = params;
problem = trinet_stack.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name', 'saigup_trinet_stacked_gdz', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_stack_tuned, p, h] = optimizeNetworkModel(trinet_stack, p0, 'maxIt', 40);
