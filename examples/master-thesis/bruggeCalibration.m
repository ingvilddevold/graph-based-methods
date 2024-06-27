mrstModule add ad-core ad-props ad-blackoil
mrstModule add deckformat
mrstModule add test-suite
mrstModule add ensemble 
mrstModule add optimization
mrstModule add mrst-gui
mrstModule add graph-based-methods
mrstModule add upscaling
mrstModule add coarsegrid
mrstModule add upr
mrstModule add diagnostics

gravity reset on

%% Set up Brugge model

[problem0, deck] = setupBrugge(); % sets up model from ECLIPSE deck
    
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
% shorter well names, BR-I-5 -> I5
scheduleTrue = shortenWellNames(scheduleTrue, @(name) erase(name(4:end), '-'));

problem0.SimulatorSetup.schedule = scheduleTrue;

%clearPackedSimulatorOutput(problem0);
simulatePackedProblem(problem0);
[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem0);

%% Parameter configuration and objective function

config = {...
     %name           include  scaling    boxlims lumping subset  relLims  mapto 
    'porevolume',       1,   'linear',       [],    [],    [], [0.01 100] 'node'
    'conntrans',        1,      'log',       [],    [],    [], [0.01 100] 'well'
    'transmissibility', 1,      'log',       [],    [],    [], [0.01 100] 'edge'
    'swl',              1,   'linear',  [0, .4],    [],    [],      []    'node'   % relperm scaling
    'swcr',             1,   'linear',  [0, .4],    [],    [],      []    'node'   % .
    'swu',              1,   'linear',  [.8, 1],    [],    [],      []    'node'   % .
    'sowcr',            1,   'linear',  [0, .4],    [],    [],      []    'node'   % .
    'krw',              1,   'linear',  [.2, 2],    [],    [],      []    'node'   % .
    'kro',              1,   'linear',  [.2, 2],    [],    [],      []    'node'   % .
    'sw',               1,   'linear',    [0 1],    [],    [],      [],   'node'}; % initial water saturation

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsRef);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesRef, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% CGNet2: Non-uniform partition-based network model
% a priori flow adaptivity using time-of-flight

[cgnet2, tof] = setupBruggeNonuniformPartition(problem0, statesRef{end});

% some saturations are slightly above 1 for some reason, which causes
% trouble with the box limits. 
% enforce a maximum of 1 manually:
for i = 1:cgnet2.network.numNodes()
    if cgnet2.state0.s(i,1) > 1
        cgnet2.state0.s(i,1) = 1;
    end
end

cgnet2.params = setupParameters(cgnet2, config);
cgnet2.model.fluid = rmfield(cgnet2.model.fluid, 'pcOW'); % remove capillary pressure

samples = struct('problem', {{cgnet2.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_cgnet2 = OptimizationProblem(samples, ...
        'parameters',       cgnet2.params, ...
        'name',     'master/brugge/cgnet2', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      false, ...
        'solverFunOptions',     {'scalarObjective', false});

[cgnet2_tuned, p_cgnet2_tuned, h_cgnet2] = optimizeNetworkModel(cgnet2, p_cgnet2, 'maxIt', 20);
h_cgnet2.name = 'CGNet2 upscaled';

%% Run CGNet2 with poor initial guess 

cgnet2dd = cgnet2;

% Discard the initial state info, run fully data-driven. (comparable to TriNet)
cgnet2dd.state0.s(:,1) = repmat(mean(problem0.SimulatorSetup.state0.s(:,1)), 1, cgnet2dd.network.numNodes);
cgnet2dd.state0.s(:,2) = repmat(mean(problem0.SimulatorSetup.state0.s(:,2)), 1, cgnet2dd.network.numNodes);
cgnet2dd.state0.pressure = repmat(mean(problem0.SimulatorSetup.state0.pressure), cgnet2dd.network.numNodes, 1);

cgnet2dd.params = setupParameters(cgnet2dd, config);

samples = struct('problem', {{cgnet2dd.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_cgnet2dd = OptimizationProblem(samples, ...
        'parameters',       cgnet2dd.params, ...
        'name',     'master/brugge/cgnet2_dd', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      false, ...
        'solverFunOptions',     {'scalarObjective', false});

[cgnet2dd_tuned, p_cgnet2dd_tuned, h_cgnet2dd] = optimizeNetworkModel(cgnet2dd, p_cgnet2dd, 'maxIt', 20);
h_cgnet2dd.name = 'CGNet2 poor initial';

%% CGNet1: Uniform model of comparable granularity

blockIx = partitionUI(modelTrue.G, [13, 8, 1]);
% % wells
mx = max(blockIx);
for k = 1:numel(WTrue)
    blockIx(WTrue(k).cells) = mx +k;
end
blockIx = processPartition(modelTrue.G, blockIx);
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

cgnet1 = PartitionNet(problem.SimulatorSetup.model, ...
                    problem.SimulatorSetup.schedule, ...
                    problem.SimulatorSetup.state0, ...
                    blockIx, ...
                    problem0.SimulatorSetup.model, ...
                    problem0.SimulatorSetup.schedule, ...
                    problem0.SimulatorSetup.state0);

% some saturations are slightly above 1 for some reason, which causes
% trouble with the box limits. 
% enforce a maximum of 1 manually:
for i = 1:cgnet1.network.numNodes()
    if cgnet1.state0.s(i,1) > 1
        cgnet1.state0.s(i,1) = 1;
    end
end

cgnet1.params = setupParameters(cgnet1, config);
cgnet1.model.fluid = rmfield(cgnet1.model.fluid, 'pcOW'); % remove capillary pressure

samples = struct('problem', {{cgnet1.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_cgnet1 = OptimizationProblem(samples, ...
        'parameters',  cgnet1.params, ...
        'name',   'master/brugge/cgnet1', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      false, ...
        'solverFunOptions',     {'scalarObjective', false});

[cgnet1_tuned, p_cgnet1_tuned, h_cgnet1] = optimizeNetworkModel(cgnet1, p_cgnet1, 'maxIt', 20);
h_cgnet1.name = 'CGNet1 upscaled';

%% CGNet1 with poor initial guess
cgnet1dd = cgnet1;

% Discard the initial state info, run fully data-driven. (comparable to TriNet)
cgnet1dd.state0.s(:,1) = repmat(mean(problem0.SimulatorSetup.state0.s(:,1)), 1, cgnet1dd.network.numNodes);
cgnet1dd.state0.s(:,2) = repmat(mean(problem0.SimulatorSetup.state0.s(:,2)), 1, cgnet1dd.network.numNodes);
cgnet1dd.state0.pressure = repmat(mean(problem0.SimulatorSetup.state0.pressure), cgnet1dd.network.numNodes,1);


cgnet1dd.params = setupParameters(cgnet1dd, config);

samples = struct('problem', {{cgnet1dd.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_cgnet1dd = OptimizationProblem(samples, ...
        'parameters',  cgnet1dd.params, ...
        'name',   'master/brugge/cgnet1_dd', ...
        'objective',             objFun, ...
        'setupType',       'simulation', ...
        'verboseSimulation',      true, ...
        'solverFunOptions',     {'scalarObjective', false});

[cgnet1dd_tuned, p_cgnet1dd_tuned, h_cgnet1dd] = optimizeNetworkModel(cgnet1dd, p_cgnet1dd, 'maxIt', 20);
h_cgnet1dd.name = 'CGNet1 poor initial';

%% Triangulation-based models
% To construct TriNets, we need a 2D model to start from.

% this takes a while
cProblem.SimulatorSetup = coarsenSetup(problem0.SimulatorSetup, [139 48 1]);

pts = cProblem.SimulatorSetup.model.fluid.krPts;
cProblem.SimulatorSetup.model = imposeRelpermScaling(cProblem.SimulatorSetup.model, ...
    'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
    'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
    'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
cProblem.SimulatorSetup.model = cProblem.SimulatorSetup.model.setupOperators(cProblem.SimulatorSetup.model.G, cProblem.SimulatorSetup.model.rock);

for i = 1:numel(cProblem.SimulatorSetup.schedule.control)
    cProblem.SimulatorSetup.schedule.control(i).W = collapseWells(cProblem.SimulatorSetup.schedule.control(i).W);
end

%% TriNet1: Uniform triangulation-based model
% using DistMesh

trinet1 = TriNet(cProblem.SimulatorSetup.model, ...
                 cProblem.SimulatorSetup.schedule, ...
                 cProblem.SimulatorSetup.state0, ...
                 'edgeFac', 0.055, ...
                 'useNodes', true);
trinet1.model.operators.pv = trinet1.model.operators.pv * sum(modelTrue.operators.pv) / sum(trinet1.model.operators.pv);
trinet1.model.fluid = rmfield(trinet1.model.fluid, 'pcOW'); % remove capillary pressure

trinet1.params = setupParameters(trinet1, config);

samples = struct('problem', {{trinet1.getPackedSimulationProblem()}}, ...
                 'num',     1);

p1 = OptimizationProblem(samples, ...
            'parameters', trinet1.params, ...
            'name',        'master/brugge/trinet1_static', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});
%p1.reset('prompt', false);

[trinet1_tuned, p1_tuned, h_trinet1] = optimizeNetworkModel(trinet1, p1, 'maxIt', 20);
h_trinet1.name = 'TriNet1';

%% TriNet2: Non-uniform triangulation-based model
% using DistMesh, with higher resolution closer to wells.

trinet2 = TriNet(cProblem.SimulatorSetup.model, ...
                 cProblem.SimulatorSetup.schedule, ...
                 cProblem.SimulatorSetup.state0, ...
                 'adaptToWells', true, ...
                 'edgeFac', 0.05);
trinet2.model.operators.pv = trinet2.model.operators.pv * sum(modelTrue.operators.pv) / sum(trinet2.model.operators.pv);
trinet2.model.fluid = rmfield(trinet2.model.fluid, 'pcOW'); % remove capillary pressure

trinet2.params = setupParameters(trinet2, config);

samples = struct('problem', {{trinet2.getPackedSimulationProblem()}}, ...
                 'num',     1);

p2 = OptimizationProblem(samples, ...
            'parameters', trinet2.params, ...
            'name',        'master/brugge/trinet2_static', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});
%p2.reset('prompt', false);

[trinet2_tuned, p2_tuned, h_trinet2] = optimizeNetworkModel(trinet2, p2, 'maxIt', 20);
h_trinet2.name = 'TriNet2';


%% Automatic refinement

% start from coarse TriNet3 (47 nodes)
trinet3 = TriNet(cProblem.SimulatorSetup.model, ...
                 cProblem.SimulatorSetup.schedule, ...
                 cProblem.SimulatorSetup.state0, ...
                 'edgeFac', 0.1);
trinet3.model.operators.pv = trinet3.model.operators.pv * sum(modelTrue.operators.pv) / sum(trinet3.model.operators.pv);
trinet3.model.fluid = rmfield(trinet3.model.fluid, 'pcOW'); % remove capillary pressure

trinet3.params = setupParameters(trinet3, config);
trinet3.paramConfig = config;

samples = struct('problem', {{trinet3.getPackedSimulationProblem()}}, ...
                 'num',     1);

p3 = OptimizationProblem(samples, ...
            'parameters', trinet3.params, ...
            'name',        'master/brugge/trinet_autoref_new', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  true, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet3_tuned, p3_tuned, h_trinet3] = optimizeNetworkModelTopology(trinet3, p3, config, 'resChangeTolRel', 1e-2, 'maxIt', 100);
h_trinet3.name = 'TriNet3 auto.ref.';


%% Results figure

% load automatic refinement results
graphsLoc = fullfile(mrstOutputDirectory(), ... 
                    'master/brugge/trinet_autoref_new/graphs/');
g1 = load(fullfile(graphsLoc, 'graph_1.mat')).netmod;
g2 = load(fullfile(graphsLoc, 'graph_2.mat')).netmod;
g3 = load(fullfile(graphsLoc, 'graph_3.mat')).netmod;
g4 = load(fullfile(graphsLoc, 'graph_4.mat')).netmod;
h_trinet3 = load(fullfile(mrstOutputDirectory(), ...
    'master/brugge/trinet_autoref_new/history.mat')).histAll;
h_trinet3.name = 'TriNet3 auto';

fig = figure('Position', [1  1 600 300]);

ax1 = subplot(2,4,1);
plotNetworkModel(g1, ax1, 'flip', true, 'label', 'none', 'wellnamesize', 8, ...
            'nodesize', 2, 'wellsize', 3);
axis off

ax2 = subplot(2,4,2);
plotNetworkModel(g2, ax2, 'flip', true, 'label', 'none', 'nodesize', 2, 'wellsize', 3);
axis off

ax3 = subplot(2,4,3);
plotNetworkModel(g3, ax3, 'flip', true, 'label', 'none', 'nodesize', 2, 'wellsize', 3);
axis off

ax4 = subplot(2,4,4);
plotNetworkModel(g4, ax4, 'flip', true, 'label', 'none', 'nodesize', 2, 'wellsize', 3);
axis off

ax = subplot(2,4,5:8);
h_cgnet1.its = h_trinet3.its; %swap order to get the right colors, first needs its.
H = {h_cgnet1, h_cgnet2, h_cgnet1_dd, h_cgnet2_dd, h_trinet1, h_trinet2, h_trinet3};
plotMismatch(H, ax, 'verticalLines', true)
ax.Position = [0.05 0.11 0.9 0.5];

ax1.Position = [0.05 0.58 0.25 0.4];
ax2.Position = [0.275 0.58 0.25 0.4];
ax3.Position = [0.5 0.58 0.25 0.4];
ax4.Position = [0.725 0.58 0.25 0.4];
