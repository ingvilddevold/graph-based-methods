%% CGNet2 control optimization
% Optimize well controls by maximizing the net-present-value (NPV), 
% using CGNet2 as our proxy model. 
%
% To run this script, you should have the calibrated CGNet2, cgnet2_tuned,
% and the original Brugge problem0 in the workspace. 
% (See bruggeCalibration.m)

problem = cgnet2_tuned.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem)
simulatePackedProblem(problem);
[ws, states] = getPackedSimulatorOutput(problem);

% Parameters for net-present-value objective function
objectiveOpts  = {'OilPrice',           50/stb, ...
                  'WaterInjectionCost'   3/stb, ...
                  'WaterProductionCost', 3/stb, ...
                  'DiscountFactor',      0.1};

% Net-present-value objective function
objFun = @(model, states, schedule, varargin)NPVOW(model, states, schedule, varargin{:}, objectiveOpts{:});

npv0 = objFun(problem.SimulatorSetup.model, states, problem.SimulatorSetup.schedule);
npv0 = sum(vertcat(npv0{:})); % intitial npv

W    = problem.SimulatorSetup.schedule.control(1).W;

bnds = processBounds(W, 'rate(inj)', [10 1000]/day, ...       % target rate bounds
                        'lrat(prod)',[-500 -10]/day, ...      % target rate bounds
                        'bhp(inj)', [160 180]*barsa, ...      % upper bhp limit bounds
                        'bhp(prod)', [50 120]*barsa);         % target bhp bounds 
maps = setupSimulationControlMappings(problem.SimulatorSetup.schedule, bnds);       
objectiveScaling = npv0;
objStruct = struct('function', objFun, ...
                   'scaling', objectiveScaling);
samples = struct('problem', {{problem}}, 'num', 1);

p = OptimizationProblem(samples, ...
                        'name',     'master/brugge/cgnet2_controlopt', ...                           
                        'objective',        objStruct,  ...                
                        'maps',                  maps,  ...
                        'setupType',     'simulation',  ...
                        'verboseSimulation',     true);      
%p.reset('prompt', false);
[us, h] = p.maximizeObjective(problem, 'objChangeTol', 1e-8, 'gradTol', 1e-5, 'maxIt', 15);
problemOpt = p.updateProblemFun(problem, us); 

problemOpt.BaseName = 'brugge_cgnet_controlopt';
problemOpt.Name = 'brugge_cgnet_controlopt';
clearPackedSimulatorOutput(problemOpt)
simulatePackedProblem(problemOpt)
[wsOpt, statesOpt] = getPackedSimulatorOutput(problemOpt);

%% Fine-scale NPV
% simulate the fine-scale model with the optimized controls from CGNet2
problem0Opt = problem0;
problem0Opt = p.updateProblemFun(problem0Opt, us);
problem0Opt.Name = 'brugge_fine_controlopt';
clearPackedSimulatorOutput(problem0Opt);
simulatePackedProblem(problem0Opt)
[ws0Opt, states0Opt] = getPackedSimulatorOutput(problem0Opt);

% Fine-scale NPV with optimized well controls
npvF = objFun(problem0Opt.SimulatorSetup.model, states0Opt, problem0Opt.SimulatorSetup.schedule);
npvF = sum(vertcat(npvF{:}));

% Fine-scale NPV with original well controls
clearPackedSimulatorOutput(problem0)
simulatePackedProblem(problem0)
[ws0, states0] = getPackedSimulatorOutput(problem0);
npvF0 = objFun(problem0.SimulatorSetup.model, states0, problem0.SimulatorSetup.schedule);
npvF0 = sum(vertcat(npvF0{:}));
