mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add network-models

%% Set up SAIGUP case

saigup = TestCase('saigup_wo');
saigup.state0.s = saigup.state0.s.*0 + [0,1];

problem = saigup.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);

%% Coarse model
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [10 10 1]);

%% Plot fine and coarse grid
figure
plotGrid(problem.SimulatorSetup.model.G);
plotGrid(cProblem.SimulatorSetup.model.G,'FaceColor','r','FaceAlpha',.05,'LineWidth',2);
plotWell(problem.SimulatorSetup.model.G,problem.SimulatorSetup.schedule.control(1).W,'FontSize',10);
view(2), axis tight
axis equal

%% Setup CGNet model
cgnet = NetworkModel(cProblem.SimulatorSetup.model, ...
                     cProblem.SimulatorSetup.schedule, ...
                     cProblem.SimulatorSetup.state0);

%% Set up parameters and objective

problem = cgnet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...   
                 'num',     1);

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
        'uniformLimits', false, 'mapTo', config{k, 8});
end

% Define objective function
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

cgnet.params = params;

%% Tune parameters

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',     'saigup_cgnet', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[cgnet, p, h] = optimizeNetworkModel(cgnet, p0, 'maxIt', 40);
