mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr

%% Set up Egg model
trainEx   = TestCase('egg_wo');

problem = trainEx.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);
%test.plot(statesFine)
%plotWellSols(wellSolsFine);

%% Need a 2D model to construct graph. Each well should only have one cell

cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [25 25 1]);

pts = cProblem.SimulatorSetup.model.fluid.krPts;
cProblem.SimulatorSetup.model = imposeRelpermScaling(cProblem.SimulatorSetup.model, ...
    'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
    'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
    'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
cProblem.SimulatorSetup.model = cProblem.SimulatorSetup.model.setupOperators(cProblem.SimulatorSetup.model.G, cProblem.SimulatorSetup.model.rock);

%% Set up TriNet
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0, ...
                'useDistmesh', true, ...
                'edgeFac', 0.2, ...
                'thres', 10);

% shorter well names
for i = 1:numel(trinet.schedule.control(1).W)
    for j = 1:numel(trinet.schedule.control)
        name = trinet.schedule.control(j).W(i).name;
        name = strcat(name(1),name(end));
        trinet.schedule.control(j).W(i).name = name;
    end
end

%% Set up parameters and objective function

problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits  mapTo
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
trinet.params = params;

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

p0 = OptimizationProblem(samples, ...
            'parameters',        params, ...
            'name',        'egg_trinet', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

%% 
[trinet_tuned, p, h] = optimizeNetworkModel(trinet, p0);

