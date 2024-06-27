mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble optimization
mrstModule add mrst-gui
mrstModule add upr
mrstModule add graph-based-methods

%% Set up demo five-well case with heterogeneous permeability and porosity.

demo = TestCase('diagnostics_2d_wo', 'lognormal', true, 'barriers', false, 'pvi', 1.5);
demo.name = 'diagnostics_2d_wo_heterogeneous'; % to separate from homogeneous case

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


%% Triangle refinement of TriNet1

trinet1 = TriNet(problem0.SimulatorSetup.model, ...
                 problem0.SimulatorSetup.schedule, ...
                 problem0.SimulatorSetup.state0, ...
                 'useDistmesh', false);

trinet1.params = setupParameters(trinet1, config);
trinet1.paramConfig = config; % needed to construct new params after refinement

problem = trinet1.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_trinet1 = OptimizationProblem(samples, ...
            'parameters', trinet1.params, ...
            'name', 'master/demo_graphopt_het_trinet1', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet1_tuned, p_trinet1_tuned, h_trinet1] = optimizeNetworkModelTopology(trinet1, p_trinet1, config, ...
                'maxFrac', 0.2, ...
                'maxIt', 200, 'outerMaxIt', 6, 'totalMaxIt', 1000, 'outerObjTol', 1e-4, 'resTolAbs', 1e-9, ...
                'plotEvolution', true, 'name', "TriNet", 'saveGraphs', true, ...
                'saveWellSols', true, 'wellSolsFine', wellSolsFine, ...
                'resChangeTolRel', 0.01);

%% Tuning finer model TriNet2 for comparison

% Construct by refining all triangles of TriNet1
trinet2 = trinet1;
trinet2.paramConfig = config;
trinet2 = trinet2.refineAllTriangles();

% Make sure to conserve the total pore volume by rescaling
trinet2.model.operators.pv = trinet2.model.operators.pv ...
                               * sum(trinet1.model.operators.pv) ... % total before
                               / sum(trinet2.model.operators.pv);    % total after

trinet2.params = setupParameters(trinet2, config);


problem = trinet2.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p_trinet2 = OptimizationProblem(samples, ...
            'parameters',   trinet2.params, ...
            'name', 'master/demo_het_static_trinet2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[trinet2_tuned, p2_tuned, h_trinet2] = optimizeNetworkModel(trinet2, p_trinet2, 'maxIt', 40);
h_trinet2.name = 'TriNet2';

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
            'name', 'master/demo_het_static_cgnet1', ...
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
            'name', 'master/demo_het_static_cgnet2', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});

[cgnet2_tuned, p_cgnet2_tuned, h_cgnet2] = optimizeNetworkModel(cgnet2, p_cgnet2, 'maxIt', 40);
h_cgnet2.name = 'CGNet2';

%% Generate result figure
% plot mismatch for each case, and the graphs during refinement. 

% load automatic refinement results
graphsLoc = fullfile(mrstOutputDirectory(), ... 
                    'master/demo_graphopt_het_trinet1/graphs/');
g1 = load(fullfile(graphsLoc, 'graph_1.mat')).netmod;
g2 = load(fullfile(graphsLoc, 'graph_2.mat')).netmod;
g3 = load(fullfile(graphsLoc, 'graph_3.mat')).netmod;
g4 = load(fullfile(graphsLoc, 'graph_4.mat')).netmod;
g5 = load(fullfile(graphsLoc, 'graph_5.mat')).netmod;
g6 = load(fullfile(graphsLoc, 'graph_6.mat')).netmod;

h_trinet1 = load(fullfile(mrstOutputDirectory(), ...
    'master/demo_graphopt_het_trinet1/history.mat')).histAll;
h_trinet1.name = 'TriNet1 auto';

% Plot
figure('Position', [1  1 600 300])

ax1 = subplot(2,6,1);
plotNetworkModel(g1, ax1, 'flip', true, 'label', 'wells', 'wellnamesize', 12, ...
            'nodesize', 4, 'wellsize', 4);
axis off

ax2 = subplot(2,6,2);
plotNetworkModel(g2, ax2, 'flip', true, 'label', 'none', 'nodesize', 4, 'wellsize', 4);
axis off

ax3 = subplot(2,6,3);
plotNetworkModel(g3, ax3, 'flip', true, 'label', 'none', 'nodesize', 4, 'wellsize', 4);
axis off

ax4 = subplot(2,6,4);
plotNetworkModel(g4, ax4, 'flip', true, 'label', 'none', 'nodesize', 4, 'wellsize', 4);
axis off

ax5 = subplot(2,6,5);
plotNetworkModel(g5, ax5, 'flip', true, 'label', 'none', 'nodesize', 3, 'wellsize', 3);
axis off

ax6 = subplot(2,6,6);
plotNetworkModel(g6, ax6, 'flip', true, 'label', 'none', 'nodesize', 2, 'wellsize', 2);
axis off

ax = subplot(2,6,7:12);
plotMismatch({h_trinet1, h_trinet2, h_cgnet1, h_cgnet2}, gca, 'verticalLines', true);

