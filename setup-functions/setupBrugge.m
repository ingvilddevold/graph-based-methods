function [problem, deck] = setupBrugge()
mfl = mfilename('fullpath');
dd = fileparts(mfl);
deck = readEclipseDeck(fullfile(dd, 'brugge_new','BRUGGE60K_FY-SF-KM-1-1.DATA'));
deck = convertDeckUnits(deck);
G = computeGeometry(initEclipseGrid(deck));
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

rock.poro(rock.poro==0) = 1e-3;
rock.ntg(rock.ntg==0) = 1;

model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
W = processWells(G, rock, deck.SCHEDULE.control(end)); 
for k =1:numel(W)
    W(k).vfp_index = 0;
    W(k).compi = W(k).compi(1:2);
end

dT = rampupTimesteps(10*year, 30*day);
schedule = simpleSchedule(dT, 'W', W);

regions = getInitializationRegionsDeck(model, deck);

state0 = initStateBlackOilAD(model, regions);

nonlinear = NonLinearSolver();
nonlinear.LinearSolver = selectLinearSolverAD(model);

problem = packSimulationProblem(state0, model, schedule, 'brugge_tune', 'NonLinearSolver', nonlinear);
end