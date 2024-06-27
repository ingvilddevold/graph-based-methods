% Calibration of a CGNet and a TriNet for the Egg model
% and testing on perturbed controls

mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble 
mrstModule add optimization
mrstModule add mrst-gui
mrstModule add graph-based-methods
mrstModule add upscaling
mrstModule add coarsegrid
mrstModule add upr

%% Set up Egg model

egg = TestCase('egg_wo');

% shorter well names, INJECT1 -> I1
egg.schedule = shortenWellNames(egg.schedule, @(name) strcat(name(1),name(end)));

% run fine-scale simulation
problem0 = egg.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0);

[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem0);
%egg.plot(statesFine)
%plotWellSols(wellSolsFine);


%% Parameter configuration and objective function

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

% Define objective function (weighted sum of squares)
weighting = objectiveWeighting(wellSolsFine);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% CGNet: Uniform partition-based network model

%--------------------
% Coarse partition
%--------------------
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

%--------------------
% Set up CGNet
%--------------------
cgnet = PartitionNet(problem.SimulatorSetup.model, ...
                    problem.SimulatorSetup.schedule, ...
                    problem.SimulatorSetup.state0, ...
                    blockIx, ...
                    problem0.SimulatorSetup.model, ...
                    problem0.SimulatorSetup.schedule, ...
                    problem0.SimulatorSetup.state0);
cgnet.name = 'egg_wo_cgnet';

%--------------------
% Calibration
%--------------------
cgnet.params = setupParameters(cgnet, config);

samples = struct('problem', {{cgnet.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_cgnet= OptimizationProblem(samples, ...
            'parameters',   cgnet.params, ...
            'name',         'master/egg/cgnet_6x6x1', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[cgnet_tuned, p_cgnet_tuned, h_cgnet] = optimizeNetworkModel(cgnet, p_cgnet, 'maxIt', 20);
h_cgnet.name = "CGNet";

%% TriNet

% Need a 2D model to construct TriNet. Each well should only have one cell
% (this takes a while)
cProblem.SimulatorSetup = coarsenSetup(problem0.SimulatorSetup, [60 60 1]);

pts = cProblem.SimulatorSetup.model.fluid.krPts;
cProblem.SimulatorSetup.model = imposeRelpermScaling(cProblem.SimulatorSetup.model, ...
    'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
    'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
    'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
cProblem.SimulatorSetup.model = cProblem.SimulatorSetup.model.setupOperators(cProblem.SimulatorSetup.model.G, cProblem.SimulatorSetup.model.rock);

% Uniform TriNet

trinet1 = TriNet(cProblem.SimulatorSetup.model, ...
                 cProblem.SimulatorSetup.schedule, ...
                 cProblem.SimulatorSetup.state0, ...
                 'useDistmesh', true, ...
                 'edgeFac', 0.13, ...
                 'thres', 10);

trinet1.model.operators.pv = trinet1.model.operators.pv * sum(egg.model.operators.pv) / sum(trinet1.model.operators.pv);

trinet1.params = setupParameters(trinet1, config);

samples = struct('problem', {{trinet1.getPackedSimulationProblem()}}, ...
                 'num',     1);

p_trinet1 = OptimizationProblem(samples, ...
            'parameters', trinet1.params, ...
            'name',    'master/egg/trinet1', ...
            'objective',         objFun, ...
            'setupType',   'simulation', ...
            'verboseSimulation',  false, ...
            'solverFunOptions', {'scalarObjective', false});
%p_trinet1.reset('prompt', false);

[trinet1_tuned, p1_tuned, h_trinet1] = optimizeNetworkModel(trinet1, p_trinet1, 'maxIt', 20);
h_trinet1.name = "TriNet";

%% Try different perturbation levels and compute objective function

ratePerturbations = [0.05 0.1 0.15 0.25];
bhpPerturbations = [5, 7.5, 10, 15];
objectivesCGNet  = [];
objectivesTriNet = [];

for i = 1:numel(ratePerturbations)
    ratePert = ratePerturbations(i);
    bhpPert = bhpPerturbations(i);
    
    % Get perturbed fine-scale model
    eggPert = makeRandomTraining(egg, ratePert, @(x,y) y - bhpPert*(x-.2)*barsa, false);
    eggPert.name = strcat(eggPert.name, num2str(i));
    % and simulate
    problem0Pert = eggPert.getPackedSimulationProblem();
    simulatePackedProblem(problem0Pert);
    [wellSolsFinePert, statesFinePert] = getPackedSimulatorOutput(problem0Pert);

    % Get perturbed CGNet
    cgnetPert = cgnet_tuned;
    cgnetPert.name = strcat(cgnetPert.name, num2str(i));
    for j = 1:numel(eggPert.schedule.control)
        for w = 1:cgnetPert.numWells
            cgnetPert.schedule.control(j).W(w) = cgnetPert.schedule.control(1).W(w);
            cgnetPert.schedule.control(j).W(w).val = eggPert.schedule.control(j).W(w).val;
        end
    end
    cgnetPert.schedule.step = eggPert.schedule.step;

    % Run single iteration of optimization to get objective value
    objectivesCGNet(i) = evaluateSingleObjective(cgnetPert, statesFinePert, weighting);

    % Get perturbed TriNet
    trinetPert = trinet1_tuned;
    trinetPert.name = strcat(trinetPert.name, num2str(i));
    for j = 1:numel(eggPert.schedule.control)
        for w = 1:trinetPert.numWells
            trinetPert.schedule.control(j).W(w) = trinetPert.schedule.control(1).W(w);
            trinetPert.schedule.control(j).W(w).val = eggPert.schedule.control(j).W(w).val;
        end
    end
    trinetPert.schedule.step = eggPert.schedule.step;

    % Run single iteration of optimization to get objective value
    objectivesTriNet(i) = evaluateSingleObjective(trinetPert, statesFinePert, weighting);
end

%% Test on perturbed schedule and plot well curves

% Perturbed fine-scale simulation.
level = 1;
egg_pert = makeRandomTraining(egg, ...
                            ratePerturbations(level), ...
                            @(x,y) y - bhpPerturbations(level)*(x-.2)*barsa, ...
                            false);
egg_pert.name = strcat('egg_wo_pert_level', num2str(level));
problem0_pert = egg_pert.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem0_pert);
simulatePackedProblem(problem0_pert);
[wellSolsFine_pert, statesFine_pert] = getPackedSimulatorOutput(problem0_pert);

% Get tuned CGNet
cgnetPert = cgnet_tuned;
cgnetPert.name = 'egg_wo_pert';
% and switch to perturbed controls
for i = 1:numel(egg_pert.schedule.control)
    for w = 1:cgnetPert.numWells
        cgnetPert.schedule.control(i).W(w) = cgnetPert.schedule.control(1).W(w);
        cgnetPert.schedule.control(i).W(w).val = egg_pert.schedule.control(i).W(w).val;
    end
end
cgnetPert.schedule.step = egg_pert.schedule.step;

problem_pert = cgnetPert.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem_pert)
simulatePackedProblem(problem_pert);
[wellSolsCoarse_pert, statesCoarse_pert] = getPackedSimulatorOutput(problem_pert);

% Get tuned TriNet
trinetPert = trinet1_tuned;
trinetPert.name = 'egg_wo_pert';
% and switch to perturbed controls
for j = 1:numel(egg_pert.schedule.control)
    for w = 1:trinetPert.numWells
        trinetPert.schedule.control(j).W(w) = trinetPert.schedule.control(1).W(w);
        trinetPert.schedule.control(j).W(w).val = egg_pert.schedule.control(j).W(w).val;
    end
end
trinetPert.schedule.step = egg_pert.schedule.step;

problem_pert = trinetPert.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem_pert)
simulatePackedProblem(problem_pert);
[wellSolsTrinet_pert, statesTrinet_pert] = getPackedSimulatorOutput(problem_pert);

% Compare well responses in a plot
plotWellSols({wellSolsFine_pert, wellSolsCoarse_pert, wellSolsTrinet_pert}, ...
    egg.schedule.step.val, ...
    'datasetnames', {'Reference', 'CGNet', 'TriNet'}, ...
    'linestyles', {'-', '--', ':'}, ...
    'field', 'qOs', ...
    'timescale', 'years');
