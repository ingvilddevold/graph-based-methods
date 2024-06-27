function weighting = objectiveWeighting(wellSolsFine)
    % Calculate weights for the mismatch function 
    % based on fine-scale solution.
    
    qWs = getWellOutput(wellSolsFine, 'qWs');
    qOs = getWellOutput(wellSolsFine, 'qOs');
    bhp = getWellOutput(wellSolsFine, 'bhp');
    
    wqWs = 1/max(abs(qWs(:)));
    wqOs = 1/max(abs(qOs(:)));
    wbhp = 1/(max(bhp(:)) - min(bhp(:)));
    
    weighting =  {'WaterRateWeight',  wqWs, ...
                  'OilRateWeight',    wqOs, ...
                  'BHPWeight',        wbhp};
end