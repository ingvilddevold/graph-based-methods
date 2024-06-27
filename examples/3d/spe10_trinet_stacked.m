mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add spe10
mrstModule add upr
mrstModule add network-models

%% Set up SPE10 case

spe10 = TestCase('spe10_wo', 'layers', 10);

problem = spe10.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);

%% Need a 2D model to construct graph. Each well should only have one cell
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [3 3 1]);

%% Set up 2D TriNet model
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0);

%% Set up parameters and objective function

problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

config = {
    ...%name      include   scaling  boxlims lumping subset relativeLimits  mapTo
    'porevolume',       1, 'linear',   [],     [],    [],     [0.01 10],   'node' 
    'conntrans',        1, 'log',      [],     [],    [],     [0.01 100],  'well'
    'transmissibility', 1, 'log'       [],     [],    [],     [0.01 100],  'edge'}; 

params = setupParameters(trinet, config);

% Define objective function
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

trinet.params = params;

%% Create model 
top = min(spe10.model.G.cells.centroids(:,3));
bottom = max(spe10.model.G.cells.centroids(:,3));
trinet_stack = stackNetworkModel(trinet, [top, bottom+10], params);

%% Static tuning

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name', 'spe10_trinet_static', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_stack_tuned, p, h] = optimizeNetworkModel(trinet_stack, p0, 'maxIt', 40);

