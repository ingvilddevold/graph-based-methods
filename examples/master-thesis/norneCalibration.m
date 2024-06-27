% Calibrate network models for Norne water-oil model
% The reservoir is initially filled with oil (Sw=0, So=1), 
% and we simulate immiscible water-flooding.
%
% We calibrate pore volumes, transmissibilities and well indices for: 
%   CGNet1 : partition-based model from 9x9x1 partition
%   CGNet2 : partition-based model from 6x6x2 partition
%   CGNet3 : partition-based model from 9x9x2 partition
%   CGNet4 : partition-based model from 6x6x1 partition

mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add optimization
mrstModule add mrst-gui
mrstModule add coarsegrid
mrstModule add graph-based-methods

%% Set up fine-scale Norne model and simulate
norne = TestCase('norne_simple_wo');

problem0 = norne.getPackedSimulationProblem();
%clearPackedSimulatorOutput();
simulatePackedProblem(problem0);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem0);

G = problem0.SimulatorSetup.model.G;

%% Visualize the model
figure
plotGrid(norne.model.G, 'EdgeAlpha', 0.1, 'FaceAlpha', 1, 'FaceColor', [0.90,0.90,0.90])
plotWell(norne.model.G, norne.schedule.control(1).W, 'color', 'k')
set(gca,'dataasp',[1 1 0.2]);
axis tight off, view([84, 56]); zoom(1), colormap(jet), axis tight

%% Parameter configuration and objective function

config = {...
     %name           include  scaling    boxlims lumping subset  relLims  mapto 
    'porevolume',       1,   'linear',       [],    [],    [], [0.01 100] 'node'
    'conntrans',        1,      'log',       [],    [],    [], [0.01 100] 'well'
    'transmissibility', 1,      'log',       [],    [],    [], [0.01 100] 'edge'};

% Define objective function
weighting = objectiveWeighting(wellSolsRef);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesRef, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% CGNet1
% Constructed from a 9x9x1 coarse partition

partition1 = partitionUI(G, [9, 9, 1]);

% wells
mx = max(partition1);
WTrue = norne.schedule.control(1).W;
for k = 1:numel(WTrue)
    partition1(WTrue(k).cells) = mx +k;
end
partition1 = processPartition(G, partition1);
partition1 = compressPartition(partition1);

modelC = upscaleModelTPFA(norne.model, partition1);
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

cgnet1 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     partition1, ...
                     problem0.SimulatorSetup.model, ...
                     problem0.SimulatorSetup.schedule, ...
                     problem0.SimulatorSetup.state0);

cgnet1.params = setupParameters(cgnet1, config);

samples = struct('problem', {{cgnet1.getPackedSimulationProblem()}}, ...
                 'num',     1);

p1 = OptimizationProblem(samples, ...
            'parameters', cgnet1.params, ...
            'name',        'master/norne/cgnet1_9x9x1', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet1_tuned, p1_tuned, h_cgnet1] = optimizeNetworkModel(cgnet1, p1, 'maxIt', 10);
h_cgnet1.name = 'CGNet 9x9x1';

%% CGNet2
% Constructed from a 6x6x2 partition

partition2 = partitionUI(G, [6, 6, 2]);

% wells
mx = max(partition2);
WTrue = norne.schedule.control(1).W;
for k = 1:numel(WTrue)
    partition2(WTrue(k).cells) = mx +k;
end
partition2 = processPartition(G, partition2);
partition2 = compressPartition(partition2);

modelC = upscaleModelTPFA(norne.model, partition2);
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

cgnet2 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     partition2, ...
                     problem0.SimulatorSetup.model, ...
                     problem0.SimulatorSetup.schedule, ...
                     problem0.SimulatorSetup.state0);

cgnet2.params = setupParameters(cgnet2, config);

samples = struct('problem', {{cgnet2.getPackedSimulationProblem()}}, ...
                 'num',     1);

p2 = OptimizationProblem(samples, ...
            'parameters', cgnet2.params, ...
            'name',        'master/norne/cgnet2_6x6x2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet2_tuned, p2_tuned, h_cgnet2] = optimizeNetworkModel(cgnet2, p2, 'maxIt', 10);
h_cgnet2.name = 'CGNet 6x6x2';

%% CGNet3
% Constructed from a 9x9x2 uniform partition

partition3 = partitionUI(G, [9, 9, 2]);

% wells
mx = max(partition3);
WTrue = norne.schedule.control(1).W;
for k = 1:numel(WTrue)
    partition3(WTrue(k).cells) = mx +k;
end
partition3 = processPartition(G, partition3);
partition3 = compressPartition(partition3);

modelC = upscaleModelTPFA(norne.model, partition3);
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

cgnet3 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     partition3, ...
                     problem0.SimulatorSetup.model, ...
                     problem0.SimulatorSetup.schedule, ...
                     problem0.SimulatorSetup.state0);

cgnet3.params = setupParameters(cgnet3, config);

samples = struct('problem', {{cgnet3.getPackedSimulationProblem()}}, ...
                 'num',     1);

p3 = OptimizationProblem(samples, ...
            'parameters', cgnet3.params, ...
            'name',        'master/norne/cgnet3_9x9x2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet3_tuned, p3_tuned, h_cgnet3] = optimizeNetworkModel(cgnet3, p3, 'maxIt', 10);
h_cgnet3.name = 'CGNet 9x9x2';

%% CGNet4
% Constructed from a 6x6x1 uniform partition

partition4 = partitionUI(G, [6, 6, 1]);

% wells
mx = max(partition4);
WTrue = norne.schedule.control(1).W;
for k = 1:numel(WTrue)
    partition4(WTrue(k).cells) = mx +k;
end
partition4 = processPartition(G, partition4);
partition4 = compressPartition(partition4);

modelC = upscaleModelTPFA(norne.model, partition4);
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

cgnet4 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     partition4, ...
                     problem0.SimulatorSetup.model, ...
                     problem0.SimulatorSetup.schedule, ...
                     problem0.SimulatorSetup.state0);

cgnet4.params = setupParameters(cgnet4, config);

samples = struct('problem', {{cgnet4.getPackedSimulationProblem()}}, ...
                 'num',     1);

p4 = OptimizationProblem(samples, ...
            'parameters', cgnet4.params, ...
            'name',        'master/norne/cgnet4_6x6x1', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet4_tuned, p4_tuned, h_cgnet4] = optimizeNetworkModel(cgnet4, p4, 'maxIt', 10);
h_cgnet4.name = 'CGNet 6x6x1';

%% Compare calibration results in a plot

H = {h_cgnet1, h_cgnet2, h_cgnet3, h_cgnet4};
plotMismatch(H)

%% Plot the four partitions side-by-side

figure('Position', [1 1 1000 250])

ax1 = subplot(1,4,1);
plotCellData(G, partition1, 'EdgeAlpha', 0)
colormap(colorcube)
view(81, 71)
axis tight off

ax2 = subplot(1,4,2);
plotCellData(G, partition2, 'EdgeAlpha', 0)
colormap(colorcube)
view(81, 71)
axis tight off

ax3 = subplot(1,4,3);
plotCellData(G, partition3, 'EdgeAlpha', 0)
colormap(colorcube)
view(81, 71)
axis tight off

ax4 = subplot(1,4,4);
plotCellData(G, partition4, 'EdgeAlpha', 0)
colormap(colorcube)
view(81, 71)
axis tight off

ax1.Position = [0 -0.3 0.25 1.3];
ax2.Position = [0.25 -0.3 0.25 1.3];
ax3.Position = [0.5 -0.3 0.25 1.3];
ax4.Position = [0.75 -0.3 0.25 1.3];
