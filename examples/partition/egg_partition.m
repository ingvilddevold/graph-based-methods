mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble
mrstModule add mrst-gui
mrstModule add upr
mrstModule add coarsegrid

%% Set up Egg model
egg = TestCase('egg_wo');

% run fine-scale simulation
problem0 = egg.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);
%egg.plot(statesFine)
%plotWellSols(wellSolsFine);

%% Coarse model

G = problem0.SimulatorSetup.model.G;
blockIx = partitionUI(G, [6, 6, 1]);

% wells
mx = max(blockIx);
WTrue = egg.schedule.control(1).W;
for k = 1:numel(WTrue)
    blockIx(WTrue(k).cells) = mx +k;
end
blockIx = processPartition(G, blockIx);
blockIx = compressPartition(blockIx);

modelC = upscaleModelTPFA(egg.model, blockIx);
modelC.AutoDiffBackend = AutoDiffBackend();

pts = modelC.fluid.krPts;
scaling = {'SWL',   pts.w(1,1), 'SWCR', pts.w(1,2), 'SWU', pts.w(1,3), ...
           'SOWCR', pts.ow(1,2), 'KRW',  pts.w(1,4), 'KRO', pts.ow(1,4)};
modelC = imposeRelpermScaling(modelC, scaling{:});

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
problem.BaseName = 'egg_wo';
problem.Name = 'egg_wo';

%%
figure
explosionView(egg.model.G, blockIx);

%% Set up network model

pnet = PartitionNet(problem.SimulatorSetup.model, ...
                    problem.SimulatorSetup.schedule, ...
                    problem.SimulatorSetup.state0, ...
                    blockIx, ...
                    problem0.SimulatorSetup.model, ...
                    problem0.SimulatorSetup.schedule, ...
                    problem0.SimulatorSetup.state0);
pnet.name = 'egg_wo_pnet';

% shorter well names
for i = 1:numel(pnet.schedule.control(1).W)
    for j = 1:numel(pnet.schedule.control)
        name = pnet.schedule.control(j).W(i).name;
        name = strcat(name(1),name(end));
        pnet.schedule.control(j).W(i).name = name;
    end
end

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
    'transmissibility', 1,   'log'        [],       [],     [],        [0.01 100],  'edge'
    'swl',              1,   'linear',    [0 .2],       [],     [],        [],      'node'
    'swcr',             1,   'linear',    [0 .3],       [],     [],        [],      'node'       
    'swu',              1,   'linear',    [.7 1],       [],     [],        [],      'node'
    'sowcr',            1,   'linear',    [0 .3],       [],     [],        [],      'node'
    'krw',              1,   'linear',    [.5 2],       [],     [],        [],      'node'
    'kro',              1,   'linear',    [.5 2],       [],     [],        [],      'node'}; 

params = setupParameters(pnet, config);

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);

objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);


p0 = OptimizationProblem(samples, ...
            'parameters',       params, ...
            'name',         'egg_pnet_6x6x1', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

pnet.params = params;

%%
[pnet_tuned, p, h] = optimizeNetworkModel(pnet, p0);


%% Test on perturbed schedule
% Run a test with perturbed the controls to see if the model is 
% somewhat general.

% Perturbed fine-scale simulation.
egg_pert1 = makeRandomTraining(egg, 0.25, @(x,y) y - 5*(x-.2)*barsa, false);
egg_pert1.name = 'egg_wo_pert1';
problem0_pert1 = egg_pert1.getPackedSimulationProblem();
simulatePackedProblem(problem0_pert1);
[wellSolsFine_pert1, statesFine_pert1] = getPackedSimulatorOutput(problem0_pert1);

% Perturbed coarse-scale simulation.
pnet_pert1 = pnet_tuned;
pnet_pert1.name = 'egg_wo_pert1';

for i = 1:numel(egg_pert1.schedule.control)
    for w = 1:pnet_pert1.numWells
        pnet_pert1.schedule.control(i).W(w) = pnet_pert1.schedule.control(1).W(w);
        pnet_pert1.schedule.control(i).W(w).val = egg_pert1.schedule.control(i).W(w).val;
    end
end
pnet_pert1.schedule.step = egg_pert1.schedule.step;

problem_pert1 = pnet_pert1.getPackedSimulationProblem();
simulatePackedProblem(problem_pert1);
[wellSolsCoarse_pert1, statesCoarse_pert1] = getPackedSimulatorOutput(problem_pert1);

% Compare well responses in a plot
plotWellSols({wellSolsFine_pert1, wellSolsCoarse_pert1}, 'datasetnames', {'Fine', 'Coarse'});


%% Try different perturbation levels and compute objective function

perturbationLevels = [0.05 0.1 0.15 0.25 0.3 0.5];
objectives  = [];

for i = 1:numel(perturbationLevels)
    pert = perturbationLevels(i);
    
    % Get perturbed fine-scale model
    eggPert = makeRandomTraining(egg, pert, @(x,y) y - 5*(x-.2)*barsa, false);
    eggPert.name = strcat(eggPert.name, num2str(pert));
    % and simulate
    problem0Pert = eggPert.getPackedSimulationProblem();
    simulatePackedProblem(problem0Pert);
    [wellSolsFinePert, statesFinePert] = getPackedSimulatorOutput(problem0Pert);

    % Get perturbed coarse-scale model
    pnetPert = pnet_tuned;
    pnetPert.name = strcat(pnetPert.name, num2str(pert));
    for j = 1:numel(eggPert.schedule.control)
        for w = 1:pnetPert.numWells
            pnetPert.schedule.control(j).W(w) = pnetPert.schedule.control(1).W(w);
            pnetPert.schedule.control(j).W(w).val = eggPert.schedule.control(j).W(w).val;
        end
    end
    pnetPert.schedule.step = eggPert.schedule.step;
    % and simulate
    %problemPert = pnetPert.getPackedSimulationProblem();
    %simulatePackedProblem(problemPert);
    %[wellSolsCoarsePert, statesCoarsePert] = getPackedSimulatorOutput(problemPert);

    % Run single iteration of optimization to get objective value
    objectives(i) = evaluateSingleObjective(pnetPert, statesFinePert, weighting)
end
