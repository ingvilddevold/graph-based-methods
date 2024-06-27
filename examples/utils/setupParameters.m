function params = setupParameters(netmod, config)
    % Set up model parameters from config, matching
    % the model dimensions fo the network model netmod
    
    if isempty(config)
        warning("config is empty. No parameters constructed.")
    end
    params = [];
    for k = 1:size(config,1)
        if config{k, 2} == 0, continue, end     % include = 0
        params = addNetworkModelParameter(params, netmod, ...
            'name',    config{k,1}, 'scaling', config{k,3}, ...
            'boxLims', config{k,4}, 'lumping', config{k,5}, ...
            'subset',  config{k,6}, 'relativeLimits',config{k,7}, ...
            'uniformLimits', false, 'mapTo', config{k,8});
    end
end