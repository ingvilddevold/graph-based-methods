function W2 = collapseWells(W)
    % given a well struct W for an original 3D model, we 
    % construct a new struct with only 1 perforation, which 
    % can be used with a 2D Network Model.

    W2 = W;

    for i = 1:numel(W)
        well = W(i);

        well.cells = well.cells(1);
        well.r = well.r(1);
        well.dir = well.dir(1);
        well.rR = well.rR(1);
        well.WI = sum(well.WI);
        well.dZ = 0;
        well.cstatus = well.cstatus(1);
        %well.cell_origin = well.cell_origin(1);

        W2(i) = well;
    end
end
