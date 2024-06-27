function [pnet, T] = setupBruggeNonuniformPartition(problem0, state)
    % Set up a non-uniform partition-based model for Brugge
    % REQUIRED PARAMETERS:
    %   problem0 - Fine-scale simulation problem.
    %   state    - State for time-of-flight computation, e.g.,
    %              statesRef{end} from fine-scale simulation
    
    G = problem0.SimulatorSetup.model.G;
    rock = problem0.SimulatorSetup.model.rock;
    W = problem0.SimulatorSetup.schedule.control(1).W;

    % Compute forward and backward time-of-flight
    Tf = computeTimeOfFlight(state, G, rock, 'wells', W, 'reverse', false);
    Tb = computeTimeOfFlight(state, G, rock, 'wells', W, 'reverse', true);
    T  = Tf + Tb;

    % Start with a uniform coarse partition
    pc = partitionUI(G, [9 5 1]);
    p = compressPartition(pc);

    % Map T to coarse blocks
    Tb = nan(max(p), 1);
    for i = 1:max(p)
        Tb(i) = sum(T(p==i));
    end
    tfac = median(Tb);
    Ib = Tb < tfac; % Ib = 1 for blocks with low residence time

    % Split high-flow blocks in 4
    
    % semi-coarse partition
    pf = partitionUI(G, [18 10 1]);
    pf = compressPartition(pf);
    % shift the partition numbers to not crash with those in original 
    % coarse partition
    pf = pf + max(p);
    
    for i = 1:numel(p) % iterate fine cells
        b = p(i);      % corresponding coarse block
        if Ib(b) == 1  % corresponding block was selected
            p(i) = pf(i); % assign cell to block in semi-coarse partition
        end
    end
    p = compressPartition(p);

    % Add separate blocks for each well:
    mx = max(p);
    for k = 1:numel(W)
        p(W(k).cells) = mx +k;
    end
    p = processPartition(G, p);
    p = compressPartition(p);
    
    %------------------------------------------------
    % Set up network model
    modelC = upscaleModelTPFA(problem0.SimulatorSetup.model, p);
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
    
    pnet = PartitionNet(problem.SimulatorSetup.model, ...
                        problem.SimulatorSetup.schedule, ...
                        problem.SimulatorSetup.state0, ...
                        p, ...
                        problem0.SimulatorSetup.model, ...
                        problem0.SimulatorSetup.schedule, ...
                        problem0.SimulatorSetup.state0);
end


