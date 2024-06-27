% Calibration of a stacked TriNet model for SAIGUP, with and without 
% including gdz (gravitational effect) as a tunable parameter.
% We include sw (initial water saturation), pore volumes,
% transmissibilities and well indices in both cases.

mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite
mrstModule add ensemble optimization
mrstModule add mrst-gui
mrstModule add upr
mrstModule add graph-based-methods

%% Set up SAIGUP model

saigup = TestCase('saigup_wo');
problem = saigup.getPackedSimulationProblem();
%clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);
[wellSolsFine, statesFine] = getPackedSimulatorOutput(problem);

% each well should only have one cell, so we first coarsen the model
cProblem.SimulatorSetup = coarsenSetup(problem.SimulatorSetup, [25 25 1]);

%% Parameter config and objective function

% set box limits for gdz to g*max(dz)
bottom = max(saigup.model.G.cells.centroids(:,3));
top = min(saigup.model.G.cells.centroids(:,3));
lim = (bottom-top)*10;

config = {
    ...%name      include   scaling  boxlims lumping subset relativeLimits  mapTo
    'porevolume',       1, 'linear',   [],     [],    [],     [0.01 10],   'node' 
    'conntrans',        1, 'log',      [],     [],    [],     [0.01 100],  'well'
    'transmissibility', 1, 'log'       [],     [],    [],     [0.01 100],  'edge'
    'sw',               1, 'linear',  [0 1],   [],    [],     [],          'node'
    'gdz',              1, 'linear', [-lim lim], [],  [],     [],          'node'}; 


weighting = objectiveWeighting(wellSolsFine);
objFun = @(model, states, schedule, varargin) ...
    matchObservedOW(model, states, schedule, statesFine, varargin{:}, ...
                    weighting{:}, 'mismatchSum', false);

%% Construct a stacked TriNet

% ... by first constructing a coarse 2D TriNet
trinet = TriNet(cProblem.SimulatorSetup.model, ...
                cProblem.SimulatorSetup.schedule, ...
                cProblem.SimulatorSetup.state0);

params = setupParameters(trinet, config(1:4,:)); % excluding gdz

% ... and stacking it at the bottom and top of the reservoir

trinet_stack = stackNetworkModel(trinet, [bottom top], params);
                  
% set up parameters for trinet_stack dimensions and add gdz manually
params_stack = setupParameters(trinet_stack, config(1:4,:));

trinet_stack.model.operators.gdz = trinet_stack.model.getGravityGradient();
params_stack_gdz = addNetworkModelParameter(params_stack, trinet_stack, ...
            'name',    config{5,1}, 'scaling', config{5,3}, ...
            'boxLims', config{5,4}, 'lumping', config{5,5}, ...
            'subset',  config{5,6}, 'relativeLimits',config{5,7}, ...
            'uniformLimits', false, ...
            'belongsTo', 'model', ...
            'location', {'operators' 'gdz'}, ...
            'mapTo', 'edge');
trinet_stack.params = params_stack;

trinet_stack.size

%% Construct a finer 2D TriNet
% The stacked model has 48 nodes and 142 edges (266 tunable params excl 
% gdz), so we try to make a 2D model of comparable size.

trinet_fine = TriNet(cProblem.SimulatorSetup.model, ...
                     cProblem.SimulatorSetup.schedule, ...
                     cProblem.SimulatorSetup.state0, ...
                     'edgeFac', 0.1);

params_fine = setupParameters(trinet_fine, config(1:4,:));
trinet_fine.params = params_fine;

trinet_fine.size

%% Calibrate the stacked TriNet
% excluding gdz
problem = trinet_stack.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0_stack = OptimizationProblem(samples, ...
            'parameters', params_stack, ...
            'name', 'saigup_trinet_stack_sw', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_stack_tuned, p_stack, h_stack] = optimizeNetworkModel(trinet_stack, p0_stack, 'maxIt', 20);
h_stack.name = "TriNet stack";

%% .. and the finer 2D TriNet

problem = trinet_fine.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0_fine = OptimizationProblem(samples, ...
            'parameters', params_fine, ...
            'name',    'saigup_trinet_fine_sw', ...
            'objective',        objFun, ...
            'setupType',  'simulation', ...
            'verboseSimulation', false, ...
            'solverFunOptions',  {'scalarObjective', false});

[trinet_fine_tuned, p_fine, h_fine] = optimizeNetworkModel(trinet_fine, p0_fine, 'maxIt', 20);
h_fine.name = "TriNet fine";

%% Calibrate the stacked TriNet with gdz
%
% WARNING
% To run this, you need to modify line 13 in PhasePotentialDifference.m
% (mrst-autodiff/ad-core/statefunctions/flux/PhasePotentialDifference.m)
% 
% from
%   gp.hasGravity = norm(model.getGravityGradient(), inf) > 0;
% to
%   gp.hasGravity = norm(value(model.getGravityGradient()), inf) > 0;
%
% since model.operators.gdz will now be initialized to an ADI variable. 

trinet_stack_gdz = trinet_stack;
trinet_stack_gdz.params = params_stack_gdz;

problem = trinet_stack.getPackedSimulationProblem();
samples = struct('problem', {{problem}}, ...
                 'num',     1);

p0_stack_gdz = OptimizationProblem(samples, ...
                'parameters', params_stack_gdz, ...
                'name', 'saigup_trinet_stack_sw_gdz', ...
                'objective',        objFun, ...
                'setupType',  'simulation', ...
                'verboseSimulation', false, ...
                'solverFunOptions',  {'scalarObjective', false});

[trinet_stack_gdz_tuned, p_stack_gdz, h_stack_gdz] = optimizeNetworkModel(trinet_stack_gdz, p0_stack_gdz, 'maxIt', 20);
h_stack_gdz.name = "TriNet stack incl. gdz";

%% Generate results figure

fig = figure;
fig.Position = [0 0 650 200];

ax1 = subplot(1,6,1);
plotNetworkModel(trinet, ax1, 'label', 'none');

ax2 = subplot(1,6,2);
plotNetworkModel(trinet_fine, ax2, 'label', 'none');

ax3 = subplot(1,6,3:6);
plotMismatch({h_stack, h_fine, h_stack_gdz}, ax3);