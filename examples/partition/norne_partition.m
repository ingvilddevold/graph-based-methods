mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr
mrstModule add deckformat
mrstModule add coarsegrid

%%
norne = TestCase('norne_simple_wo');

problem0 = norne.getPackedSimulationProblem();
%clearPackedSimulatorOutput();
simulatePackedProblem(problem0);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);

% simulation not working for norne_field_bo
%  Exception thrown: The logical indices contain a true value outside of the array bounds.


%% Coarse model

G = problem0.SimulatorSetup.model.G;
blockIx = partitionUI(G, [6, 6, 1]);

% wells
mx = max(blockIx);
WTrue = norne.schedule.control(1).W;
for k = 1:numel(WTrue)
    blockIx(WTrue(k).cells) = mx +k;
end
blockIx = processPartition(G, blockIx);
blockIx = compressPartition(blockIx);

modelC = upscaleModelTPFA(norne.model, blockIx);
modelC.AutoDiffBackend = AutoDiffBackend();

% initial state/schedule
stateC0 = upscaleState(modelC, problem0.SimulatorSetup.model, ...
                               problem0.SimulatorSetup.state0);
scheduleC = upscaleSchedule(modelC, problem0.SimulatorSetup.schedule, ...
                            'wellUpscaleMethod', 'sum');

modelC.operators.T = max(modelC.operators.T, 1e-13);

setupC = struct('model', modelC, ...
                'schedule', scheduleC, ...
                'state0', stateC0);

problem.SimulatorSetup = setupC;
problem.BaseName = 'norne_simple_wo';
problem.Name = 'norne_simple_wo';

%% Plot partition
figure
explosionView(norne.model.G, blockIx);

%% Set up network model
problem0.SimulatorSetup.state0 = rmfield(problem0.SimulatorSetup.state0, 'wellSol');
problem.SimulatorSetup.state0 = rmfield(problem.SimulatorSetup.state0, 'wellSol');


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
problem = pnet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);
config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits  mapTo
    'porevolume',       1,   'linear',    [],       [],     [],        [0.01 10],   'node'
    'conntrans',        1,   'log',       [],       [],     [],        [0.01 100],  'well'
    'transmissibility', 1,   'log'        [],       [],     [],        [0.01 100],  'edge'}; 

params = setupParameters(pnet, config);

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',         'norne_pnet_7x7x1', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

pnet.params = params;

%%
[pnet_tuned, p, h] = optimizeNetworkModel(pnet, p0);
