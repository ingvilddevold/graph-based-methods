% Calibration of a coarse TriNet model for SAIGUP, with and without 
% including sw (initial water saturation) as a tunable parameter.

mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble 
mrstModule add optimization
mrstModule add mrst-gui
mrstModule add upr
mrstModule add graph-based-methods

%% Set up SAIGUP model

saigup = TestCase('saigup_wo');
problem0 = saigup.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0);
[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);

%% Construct 2D TriNet

% each well should only have one cell, so we first coarsen the model
cProblem.SimulatorSetup = coarsenSetup(problem0.SimulatorSetup, [25 25 1]);

%%
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0, ...
                'useDistmesh', true, ...
                'useNodes', true, ...
                'thres', 100, ...
                'fixBoundary', true, ...
                'edgeFac', 0.25);

%% Set up parameters and objective function

config = {
    ...%name      include   scaling  boxlims lumping subset relativeLimits  mapTo
    'porevolume',       1, 'linear',   [],     [],    [],     [0.01 10],   'node' 
    'conntrans',        1, 'log',      [],     [],    [],     [0.01 100],  'well'
    'transmissibility', 1, 'log'       [],     [],    [],     [0.01 100],  'edge'
    'sw',               1, 'linear', [0 1],    [],    [],     [],          'node'}; 

params_sw = setupParameters(trinet, config); % include sw
params = params_sw(1:3);                     % exclude sw
                                
weighting = objectiveWeighting(wellSolsFine);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% Calibrate 2D TriNet

% without saturation (sw)

problem = trinet.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',    'master/saigup_trinet_excl_sw', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_tuned, p, h] = optimizeNetworkModel(trinet, p0, 'maxIt', 20);

%% 
% with saturation

p0_sw = OptimizationProblem(samples, ...
            'parameters',       params_sw, ...
            'name',    'master/saigup_trinet_incl_sw', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_tuned_sw, p_sw, h_sw] = optimizeNetworkModel(trinet, p0_sw, 'maxIt', 20);

%% generate result figure
% showing the graph and mismatch reduction for each case, the TriNet,
% and the calibrated saturation on the Voronoi grid.

fig = figure;
fig.Position = [440 584 518 270];

ax1 = subplot(1,5,4);
plotNetworkModel(trinet, ax1, 'label', 'none');

ax2 = subplot(1,5,1:3);
h.name = "Excl. sw";
h_sw.name = "Incl. sw";
plotMismatch({h, h_sw}, ax2);

ax3 = subplot(1,5,5);
voronoiG = trinet.makeGrid();
plotCellData(voronoiG, trinet_tuned_sw.state0.s(:,2))
axis equal tight off

%% plot original oil saturation and wells
figure
plotCellData(saigup.model.G, saigup.state0.s(:,2), 'EdgeAlpha',1)
plotWell(saigup.model.G, saigup.schedule.control(1).W)
axis equal tight

%% plot calibrated oil saturation
voronoiG = trinet.makeGrid();
voronoiG.nodes.coords(:,3) = zeros(voronoiG.nodes.num, 1);
figure
plotCellData(voronoiG, trinet_tuned_sw.state0.s(:,2))
axis equal tight
