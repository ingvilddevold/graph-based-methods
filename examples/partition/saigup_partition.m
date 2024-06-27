mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr
mrstModule add coarsegrid

%% Set up SAIGUP case

saigup = TestCase('saigup_wo');
saigup.state0.s = saigup.state0.s.*0 + [0,1];

problem0 = saigup.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);

%% Coarse model

G = problem0.SimulatorSetup.model.G;
blockIx = partitionUI(G, [4, 11, 1]);
% wells
mx = max(blockIx);
WTrue = saigup.schedule.control(1).W;
for k = 1:numel(WTrue)
    blockIx(WTrue(k).cells) = mx +k;
end
blockIx = processPartition(G, blockIx);
blockIx = compressPartition(blockIx);

modelC = upscaleModelTPFA(saigup.model, blockIx);
modelC.AutoDiffBackend = AutoDiffBackend();

% initial state/schedule
stateC0 = upscaleState(modelC, problem0.SimulatorSetup.model, ...
                               problem0.SimulatorSetup.state0);
scheduleC = upscaleSchedule(modelC, problem0.SimulatorSetup.schedule, 'wellUpscaleMethod', 'sum');
modelC.operators.T = max(modelC.operators.T, 1e-13);
setupC = struct('model', modelC, 'schedule', scheduleC, 'state0', stateC0);

problem.SimulatorSetup = setupC;
problem.BaseName = 'tmp';
problem.Name = 'tmp';

%%
figure
explosionView(saigup.model.G, blockIx, 0.2)

%% Set up network model

pnet = PartitionNet(problem.SimulatorSetup.model, ...
                    problem.SimulatorSetup.schedule, ...
                    problem.SimulatorSetup.state0, ...
                    blockIx, ...
                    problem0.SimulatorSetup.model, ...
                    problem0.SimulatorSetup.schedule, ...
                    problem0.SimulatorSetup.state0);

%% Plot network model over fine model
figure
plotGrid(problem0.SimulatorSetup.model.G, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
plotNetworkModel(pnet, gcf, 'nodesize', 7, 'wellsize', 7, 'linewidth', 2, 'wellnamesize', 16)

%%
config = {...
     %name           include  scaling    boxlims lumping subset relativeLimits mapto 
    'porevolume',       1,   'linear',       [],    [],    [], [0.01 10] 'node'
    'conntrans',        1,      'log',       [],    [],    [], [0.01 100] 'well'
    'transmissibility', 1,      'log',       [],    [],    [], [0.01 100] 'edge'};

params = setupParameters(pnet, config);
pnet.params = params;

% Define objective function
weighting = objectiveWeighting(wellSolsFine);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%%
problem = pnet.getPackedSimulationProblem('saigup_wo');
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
        'parameters',            params, ...
        'name',           'saigup_pnet', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      false, ...
        'solverFunOptions',     {'scalarObjective', false});

%%
[pnet_tuned, p, h] = optimizeNetworkModel(pnet, p0, 'maxIt', 40);
