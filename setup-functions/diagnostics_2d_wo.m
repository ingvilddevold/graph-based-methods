function setup = diagnostics_2d_wo(varargin)
% Setup function for a quarter five-spot water-oil model
%
% SYNOPSIS:
%   setup = diagnostics_2d_wo('pn1', pv1, ...)
%   setup = diagnostics_2d_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a standard quarter of a five-spot case for an immisicble
%   two-phase fluid. The configurable parameters are documented below
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    % One-line description
    description = ['Two-phase oil/water version of interactiveSimple ', ...
        'example from the flow diagnostics chapter (13.5.1) of the '  , ...
        'MRST book (doi: 10.1017/9781108591416)'                      ];
    
    % Optional input arguments
    options = struct( ...
        'ncells'   , [64,64] , ... % Number of cells in x- and y-directions
        'time'     , 2*year  , ... % Total injection time
        'dt'       , 30*day  , ... % Timestep length
        'pvi'      , 1       , ... % Pore volumes injected
        'nkr'      , [2,2]   , ... % Brooks-Corey relperm exponents
        'poro'     , 0.4     , ... % Mean porosity
        'logNormal', true    , ... % Use log-normal porosity field
        'seed'     , 20221018, ... % Seed for random numbers
        'barriers' , false     ... % Include barriers on blocking flow
     ); 
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil ensemble

    % Cartesian grid covering a 1000 x 1000 m domain 
    G = computeGeometry(cartGrid(options.ncells, [5000, 2500]*meter));
    
    % Make porosity field
    if options.logNormal
        % ... from log-normal distribution with mean options.poro
        rng(options.seed); % Reproducibility
        fun = @(x) exp(-sqrt(sum(x.^2,2))/0.3);
        p = GaussianProcessND(options.ncells, fun);
        poroSTD = 0.25*options.poro;
        p = p - mean(p(:));
        p = p./std(p(:));
        poro = options.poro + p.*poroSTD;
    else
        % ... constant
        poro = repmat(options.poro, options.ncells);
    end
    
    if options.barriers
        % Add sealing barriers by setting a very low porosity
        poro(round(2*end/3):end,round(end/3)) = 1e-3;
        poro(1:round(end/3),round(2*end/3)) = 1e-3;
    end
    poro = poro(:);
    
    perm = poro.^3.*(1.5e-5)^2./(0.81*72*(1-poro).^2);
        
    rock  = makeRock(G, perm, poro); % Rock with 100 md perm and 0.4 poro
    fluid = initSimpleADIFluid('phases', 'WO'             , ... % Water and oil
                               'n'     , options.nkr      , ... % Relperm exponents
                               'mu'    , [1,1]*centi*poise, ... % Viscosity
                               'rho'   , [1,1]            );    % Density (no gravity)

    % Construct two-phase model with oil and water
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    
    % Wells
    n = 12;
    rate = options.pvi.*sum(model.operators.pv)/options.time;
    bhp  = 100*barsa;
    nx = options.ncells(1);
    W = addWell([],  G, rock, nx*(n-6)+n/2+1, ...
        'Type', 'rate', 'Comp_i', [1,0], 'name', 'I1', 'Val', rate*0.5);
    W = addWell(W, G, rock, nx*(n-6)+n/2+1+nx-n, ...
        'Type','rate',  'Comp_i', [1,0], 'name', 'I2', 'Val', rate*0.5);
    W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx-.3*nx), ...
        'Type','bhp',  'Comp_i', [0,1], 'name', 'P1', 'Val', bhp);
    W = addWell(W, G, rock, G.cells.num-(n-.5)*nx, ...
        'Type','bhp',  'Comp_i', [0,1], 'name', 'P2', 'Val', bhp);
    W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx+.3*nx), ...
        'Type','bhp',  'Comp_i', [0,1], 'name', 'P3', 'Val', bhp);
    
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(options.time, options.dt), 'W', W);
    
    % Initial state
    state0 = initResSol(G, 100*barsa, [0,1]);
    
    % Plotting
    plotOptions  = {'PlotBoxAspectRatio', [2,1,1]   , ...
                    'Size'              , [800, 400]};
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                   ...
                              'description', description , ...
                              'options'    , options     , ...
                              'state0'     , state0      , ...
                              'model'      , model       , ...
                              'schedule'   , schedule    , ...
                              'plotOptions', plotOptions);

end