function [netmod, p, h] = optimizeNetworkModel(netmod, p, varargin)
    % Optimize network model parameters.
    %
    % SYNOPSIS:
    %   [netmod, p, h] = optimizeNetworkModel(netmod, p, 'pn1', pv1, ...)
    %   
    % DESCRIPTION: 
    %
    % 
    % REQUIRED PARAMETERS:
    %   netmod - NetworkModel.
    %   p      - OptimizationProblem.
    %   
    % RETURNS:
    %   netmod  - Optimized NetworkModel.
    %   p       - Latest OptimizationProblem.
    %   h       - History struct. 
    %
    % OPTIONAL PARAMETERS:
    %   Any arguments supported by optimize and unitBoxLM, e.g. 'maxIt',
    % 

    hash = obj2hash(netmod);
    p.directory = fullfile(p.mainDirectory, hash);

    pinit = getScaledParameterVector(netmod, p.parameters);
    [popt, h] = p.optimize(pinit, 'optimizer', 'unitBoxLM', varargin{:});
    
    p.samples.problem{1} = p.updateProblemFun(p.samples.problem{1}, popt);
    netmod.model    = p.samples.problem{1}.SimulatorSetup.model;
    netmod.schedule = p.samples.problem{1}.SimulatorSetup.schedule;
    netmod.state0   = p.samples.problem{1}.SimulatorSetup.state0;
end
