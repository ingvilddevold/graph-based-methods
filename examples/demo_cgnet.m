mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui

%% Set up demo water-oil model

test = TestCase('diagnostics_2d_wo', 'lognormal', false, 'barriers', false, 'pvi', 1.5);
%test.plot('plotWells', true)

problem = test.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);
%test.plot(statesFine)
%plotWellSols(wellSols);

%% Coarse model
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [7 4]);    
cgnet = NetworkModel(cProblem.SimulatorSetup.model, ...
                     cProblem.SimulatorSetup.schedule, ...
                     cProblem.SimulatorSetup.state0);

%% Set up parameters and objective function

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
weighting =  {'WaterRateWeight',  day/3500, ...
              'OilRateWeight',    day/4500, ...
              'BHPWeight',        1/(3500*barsa)};

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);


%% Tune parameters

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name', 'diagnostics_wo_cgnet', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[cgnet, p, h] = optimizeNetworkModel(cgnet, p0, 'maxIt', 40);
