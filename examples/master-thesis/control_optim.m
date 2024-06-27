% assuming the model is in the workspace as problem0

simulatePackedProblem(problem0);
[ws0, states0] = getPackedSimulatorOutput(problem0);


% parameters for NPV
objectiveOpts  = {'OilPrice',           50/stb, ...
                  'WaterInjectionCost'   3/stb, ...
                  'WaterProductionCost', 3/stb, ...
                  'DiscountFactor',      0.1};


objFun = @(model, states, schedule, varargin)NPVOW(model, states, schedule, varargin{:}, npvWeights{:});

npv0 = objFun(problem0.SimulatorSetup.model, states0, problem0.SimulatorSetup.schedule);
npv0 = sum(vertcat(npv0{:})); % initial NPV

W      = problem0.SimulatorSetup.schedule.control(1).W;
bnds = processBounds(W, 'rate(inj)', [10 1000]/day, ...       % target rate bounds
                        'lrat(prod)',[-500 -10]/day, ...      % target rate bounds
                        'bhp(inj)', [160 180]*barsa, ...     % upper bhp limit bounds
                        'bhp(prod)', [50 120]*barsa);        % target bhp bounds 
maps = setupSimulationControlMappings(problem0.SimulatorSetup.schedule, bnds);       
objectiveScaling = npv0;
objStruct = struct('function', objFun, ...
                   'scaling', objectiveScaling);
samples = struct('problem', {{problem0}}, 'num', 1);

p = OptimizationProblem(samples, ...
                        'name',     'brugge_optimize', ...      % problem name                                
                        'objective',        objStruct,  ...                
                        'maps',                  maps,  ...
                        'setupType',     'simulation',  ...
                        'verboseSimulation',     true);      
%p.reset('prompt', false);
us = p.maximizeObjective(problem0, 'objChangeTol', 1e-8, 'gradTol', 1e-5, 'maxIt', 50);
problemOpt = p.updateProblemFun(problem0, us); 






