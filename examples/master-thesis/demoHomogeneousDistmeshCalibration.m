mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble optimization
mrstModule add mrst-gui
mrstModule add upr
mrstModule add graph-based-methods

%% Set up demo five-well case with homogeneous permeability and porosity.

demo = TestCase('diagnostics_2d_wo', 'lognormal', false, 'barriers', false, 'pvi', 1.5);
demo.name = 'diagnostics_2d_wo_homogeneous'; % to separate from heterogeneous case

problem0 = demo.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);
%demo.plot(statesFine)
%plotWellSols(wellSolsFine);

%% Parameter configuration and objective function

config = {
    ...%name      include     scaling    boxlims  lumping   subset  relativeLimits mapTo
    'porevolume',       1,   'linear',    [],       [],     [],        [0.01 10]    'node'
    'conntrans',        1,   'log',       [],       [],     [],        [0.01 100]   'well'
    'transmissibility', 1,   'log'        [],       [],     [],        [0.01 100]   'edge'}; 

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);


%% Triangle refinement of TriNet3

trinet3 = TriNet(problem0.SimulatorSetup.model, ...
                 problem0.SimulatorSetup.schedule, ...
                 problem0.SimulatorSetup.state0, ...
                 'useNodes', true, ...
                 'fixBoundary', true, ...
                 'useDistmesh', true);

trinet3.params = setupParameters(trinet3, config);
trinet3.paramConfig = config; % needed to construct new params after refinement

problem = trinet3.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_trinet3 = OptimizationProblem(samples, ...
            'parameters', trinet3.params, ...
            'name', 'master/demo_graphopt_homo_trinet3', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet3_tuned, p_trinet3_tuned, h_trinet3] = optimizeNetworkModelTopology(trinet3, p_trinet3, config, ...
                'maxFrac', 0.2, ...
                'maxIt', 200, 'outerMaxIt', 20, 'totalMaxIt', 1000, 'outerObjTol', 1e-4, 'resTolAbs', 1e-9, ...
                'plotEvolution', true, 'name', "TriNet", 'saveGraphs', true, ...
                'saveWellSols', true, 'wellSolsFine', wellSolsFine, ...
                'resChangeTolRel', 0.01);
h_trinet3.name = 'TriNet3';

%% plot TriNet3 on physical and circular layout
trinet3.plot('label', 'numbers')
trinet3.plotCircular('label', 'numbers')

%% Tuning finer model TriNet4 for comparison

trinet4 = TriNet(problem0.SimulatorSetup.model, ...
                 problem0.SimulatorSetup.schedule, ...
                 problem0.SimulatorSetup.state0, ...
                 'useDistmesh', true, ...
                 'useNodes', true, ...
                 'fixBoundary', true, ...
                 'edgeFac', 0.1);

trinet4.params = setupParameters(trinet4, config);
trinet4.paramConfig = config;

problem = trinet4.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_trinet4 = OptimizationProblem(samples, ...
            'parameters',   trinet4.params, ...
            'name', 'master/demo_homo_static_trinet4', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet4_tuned, p4_tuned, h_trinet4] = optimizeNetworkModel(trinet4, p_trinet4, 'maxIt', 40);
h_trinet4.name = 'TriNet4';

%% CGNet1, a 7x4 coarse-grid network model

cProblem1.SimulatorSetup = coarsenSetup(problem0.SimulatorSetup, [7 4]);
cgnet1 = NetworkModel(cProblem1.SimulatorSetup.model, ...
                      cProblem1.SimulatorSetup.schedule, ...
                      cProblem1.SimulatorSetup.state0);

cgnet1.params = setupParameters(cgnet1, config);

problem = cgnet1.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_cgnet1 = OptimizationProblem(samples, ...
            'parameters',   cgnet1.params, ...
            'name', 'master/demo_homo_static_cgnet1', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet1_tuned, p_cgnet1_tuned, h_cgnet1] = optimizeNetworkModel(cgnet1, p_cgnet1, 'maxIt', 40);
h_cgnet1.name = 'CGNet1';

%% CGNet2, a 13x6 coarse-grid network model

cProblem2.SimulatorSetup = coarsenSetup(problem0.SimulatorSetup, [13 6]);
cgnet2 = NetworkModel(cProblem2.SimulatorSetup.model, ...
                      cProblem2.SimulatorSetup.schedule, ...
                      cProblem2.SimulatorSetup.state0);

cgnet2.params = setupParameters(cgnet2, config);

problem = cgnet2.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_cgnet2 = OptimizationProblem(samples, ...
            'parameters',   cgnet2.params, ...
            'name', 'master/demo_homo_static_cgnet2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet2_tuned, p_cgnet2_tuned, h_cgnet2] = optimizeNetworkModel(cgnet2, p_cgnet2, 'maxIt', 40);
h_cgnet2.name = 'CGNet2';


%% Results figure
% plot mismatch for each case, and the graphs during refinement. 

% load automatic refinement results
graphsLoc = fullfile(mrstOutputDirectory(), ... 
                    'master/demo_graphopt_homo_trinet3/graphs/');
g1 = load(fullfile(graphsLoc, 'graph_1.mat')).netmod;
g2 = load(fullfile(graphsLoc, 'graph_2.mat')).netmod;
g3 = load(fullfile(graphsLoc, 'graph_3.mat')).netmod;
h_trinet3 = load(fullfile(mrstOutputDirectory(), ...
    'master/demo_graphopt_homo_trinet3/history.mat')).histAll;
h_trinet3.name = 'TriNet3 auto';

% Plot
figure('Position', [1  1 600 300])

ax1 = subplot(3,3,1);
plotNetworkModel(g1, ax1, 'flip', false, 'label', 'wells', 'wellnamesize', 12, ...
            'nodesize', 4, 'wellsize', 4);
axis off

ax2 = subplot(3,3,2);
plotNetworkModel(g2, ax2, 'flip', false, 'label', 'none', 'nodesize', 4, 'wellsize', 4);
axis off

ax3 = subplot(3,3,3);
plotNetworkModel(g3, ax3, 'flip', false, 'label', 'none', 'nodesize', 4, 'wellsize', 4);
axis off

ax = subplot(3,3,4:9);
plotMismatch({h_trinet3, h_trinet4, h_cgnet1, h_cgnet2}, gca, 'verticalLines', true);
