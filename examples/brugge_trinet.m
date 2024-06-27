%% Setup simulation problem
mrstModule add ad-core ad-blackoil deckformat diagnostics...
               mrst-gui ad-props incomp optimization...
               linearsolvers ...
               ensemble coarsegrid

gravity reset on

[problem0, deck] = setupBrugge();

% Reduce steps while debugging
%{
s = problem0.SimulatorSetup.schedule;
s.step.control = s.step.control(10);
s.step.val    = s.step.val(1:10);
problem0.SimulatorSetup.schedule = s;
%}

dt = problem0.SimulatorSetup.schedule.step.val;
nsteps = numel(dt);
perturbStep = ceil(cumsum(dt)/(year));

modelTrue    = problem0.SimulatorSetup.model;
pts = modelTrue.fluid.krPts;
scaling = {'SWL',   pts.w(1,1), 'SWCR', pts.w(1,2), 'SWU', pts.w(1,3), ...
           'SOWCR', pts.ow(1,2), 'KRW',  pts.w(1,4), 'KRO', pts.ow(1,4)};
modelTrue = imposeRelpermScaling(modelTrue, scaling{:});
problem0.SimulatorSetup.model = modelTrue;

WTrue = problem0.SimulatorSetup.schedule.control.W;
dt = problem0.SimulatorSetup.schedule.step.val;
scheduleTrue = perturbedSimpleSchedule(dt, 'W', WTrue, 'pressureFac', 0, ...
                                   'rateFac', 0, 'perturbStep', perturbStep);
problem0.SimulatorSetup.schedule = scheduleTrue;

%clearPackedSimulatorOutput(problem0);
simulatePackedProblem(problem0);
[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem0);
          
%% Coarse-scale model
G = problem0.SimulatorSetup.model.G;
blockIx = partitionUI(G, [11, 11, 1]);
% % wells
mx = max(blockIx);
for k = 1:numel(WTrue)
    blockIx(WTrue(k).cells) = mx +k;
end
blockIx = processPartition(G, blockIx);
blockIx = compressPartition(blockIx);

modelC = upscaleModelTPFA(modelTrue, blockIx);
modelC.AutoDiffBackend = AutoDiffBackend();

pts = modelC.fluid.krPts;
scaling = {'SWL',   pts.w(1,1), 'SWCR', pts.w(1,2), 'SWU', pts.w(1,3), ...
           'SOWCR', pts.ow(1,2), 'KRW',  pts.w(1,4), 'KRO', pts.ow(1,4)};
modelC = imposeRelpermScaling(modelC, scaling{:});
modelC.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

% initial state/schedule
stateC0 = upscaleState(modelC, problem0.SimulatorSetup.model, ...
                               problem0.SimulatorSetup.state0);
scheduleC = upscaleSchedule(modelC, problem0.SimulatorSetup.schedule, 'wellUpscaleMethod', 'sum');
modelC.operators.T = max(modelC.operators.T, 1e-13);
setupC = struct('model', modelC, 'schedule', scheduleC, 'state0', stateC0);
problem.SimulatorSetup = setupC;

%% Test Delaunay graph

trinet = TriNet(problem.SimulatorSetup.model, ...
                problem.SimulatorSetup.schedule, ...
                problem.SimulatorSetup.state0, ...
                'edgeFac', 0.1);

config = {...
     %name           include   scaling    boxlims  lumping subset  relativeLimits mapto 
    'porevolume',       1,    'linear',     [],    [],    [],     [0.01, 10]      'node'
    'conntrans',        1,       'log',     [],    [],    [],     [0.01 100]      'well'
    'transmissibility', 1,       'log',     [],    [],    [],     [0.01 100]      'edge'};
%{
    'swl',              1,    'linear',          [0, .4],    [],    [],      []  'node'
    'swcr',             1,    'linear',          [0, .4],    [],    [],      []  'node' 
    'swu',              1,    'linear',          [.8, 1],    [],    [],      []  'node'
    'sowcr',            1,    'linear',          [0, .4],    [],    [],      []  'node'
    'krw',              1,    'linear',        [.5, 1.5],    [],    [],      []  'node'
    'kro',              1,    'linear',        [.2, 1.5],    [],    [],      []  'node'};
%}

params = setupParameters(trinet, config);
trinet.params = params;

% Define objective function
weighting = objectiveWeighting(wellSolsRef);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesRef, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%%
problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
        'parameters',            params, ...
        'name',        'brugge_trinet', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      false, ...
        'solverFunOptions',     {'scalarObjective', false});

[trinet_tuned, p, h] = optimizeNetworkModel(trinet, p0, 'maxIt', 40);
