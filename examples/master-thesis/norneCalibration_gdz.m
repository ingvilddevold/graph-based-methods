% Calibrate network models for Norne water-oil model
% The reservoir is initially filled with oil (Sw=0, So=1)
% and we simulate immiscible water-flooding.
%
% We use the models
%   CGNet2 : partition-based model from 6*6*2 partition
%   CGNet3 : partition-based model from 9*9*1 partition
% first calibrating pore volumes, transmissibilities and well indices, and
% then including gdz (gravitational effect)

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

%% CGNet2: 3-dimensional coarse-grid network
% Constructed from a 6x6x2 partition

blockIx = partitionUI(G, [6, 6, 2]);
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

cgnet2 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     blockIx, ...
                     problem0.SimulatorSetup.model, ...
                     problem0.SimulatorSetup.schedule, ...
                     problem0.SimulatorSetup.state0);

cgnet2.params = setupParameters(cgnet2, config);

samples = struct('problem', {{cgnet2.getPackedSimulationProblem()}}, ...
                 'num',     1);

p2 = OptimizationProblem(samples, ...
            'parameters', cgnet2.params, ...
            'name',        'master/norne/cgnet2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet2_tuned, p2_tuned, h_cgnet2] = optimizeNetworkModel(cgnet2, p2, 'maxIt', 10);
h_cgnet2.name = 'CGNet 6x6x2';

%% Calibrate CGNet2 with gdz

cgnet2gdz = cgnet2;
cgnet2gdz.model.operators.gdz = cgnet2gdz.model.getGravityGradient();

% Set box lims based on top and bottom z coordinate
top = min(norne.model.G.cells.centroids(:,3));
bottom = max(norne.model.G.cells.centroids(:,3));
ulim = (bottom-top)*10;

config_gdz = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'gdz',              1,   'linear',  [-ulim, ulim],   [],     [],    [], 'edge'}; 

cgnet2gdz.params = addNetworkModelParameter(cgnet2gdz.params, cgnet2gdz, ...
        'name',    config_gdz{1,1}, 'scaling', config_gdz{1,3}, ...
        'boxLims', config_gdz{1,4}, 'lumping', config_gdz{1,5}, ...
        'subset',  config_gdz{1,6}, 'relativeLimits',config_gdz{1,7}, ...
        'uniformLimits', false, ...
        'belongsTo', 'model', ...
        'location', {'operators' 'gdz'}, ...
        'mapTo', 'edge');

p2gdz = OptimizationProblem(samples, ...
            'parameters', cgnet2gdz.params, ...
            'name',        'master/norne/cgnet2gdz', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet2gdz_tuned, p2gdz_tuned, h_cgnet2gdz] = optimizeNetworkModel(cgnet2gdz, p3, 'maxIt', 10);
h_cgnet2gdz.name = 'CGNet 6x6x2 gdz';


%% CGNet3
% Constructed from a 9x9x2 uniform partition

blockIx = partitionUI(G, [9, 9, 2]);

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

cgnet3 = PartitionNet(problem.SimulatorSetup.model, ...
                     problem.SimulatorSetup.schedule, ...
                     problem.SimulatorSetup.state0, ...
                     blockIx, ...
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

%% Calibrate CGNet3 with gdz

cgnet3gdz = cgnet3;
cgnet3gdz.model.operators.gdz = cgnet3gdz.model.getGravityGradient();

% Set box lims based on top and bottom z coordinate
top = min(norne.model.G.cells.centroids(:,3));
bottom = max(norne.model.G.cells.centroids(:,3));
ulim = (bottom-top)*10;

config_gdz = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits
    'gdz',              1,   'linear',  [-ulim, ulim],   [],     [],    [], 'edge'}; 

cgnet3gdz.params = addNetworkModelParameter(cgnet3gdz.params, cgnet3gdz, ...
        'name',    config_gdz{1,1}, 'scaling', config_gdz{1,3}, ...
        'boxLims', config_gdz{1,4}, 'lumping', config_gdz{1,5}, ...
        'subset',  config_gdz{1,6}, 'relativeLimits',config_gdz{1,7}, ...
        'uniformLimits', false, ...
        'belongsTo', 'model', ...
        'location', {'operators' 'gdz'}, ...
        'mapTo', 'edge');

p3gdz = OptimizationProblem(samples, ...
            'parameters', cgnet3gdz.params, ...
            'name',        'master/norne/cgnet3gdz', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet3gdz_tuned, p3gdz_tuned, h_cgnet3gdz] = optimizeNetworkModel(cgnet3gdz, p3gdz, 'maxIt', 10);
h_cgnet3gdz.name = 'CGNet 9x9x2 gdz';

%% Compare calibration results in mismatch plot

H = {h_cgnet2, h_cgnet3, h_cgnet2gdz, h_cgnet3gdz};
plotMismatch(H)
