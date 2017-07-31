%% quarterSpotDP - This example outlines the basic usage of the 
% two-phase dual-porosity model. We inject the first phase on 
% one corner and produce at the opposite one.
clc;
close all;

%% Include the essential modules. dual-porosity model has the functionality
% we need.
mrstModule clear
mrstModule add ad-props ad-core ad-blackoil ad-fi blackoil-sequential dual-porosity

%% Create the grid: a simple cartesian grid, 100m x 100m and 50 blocks on x
% and y directions.

% Two layer grid
x_size = 100;
y_size = 100;
z_size = 1;

Nx = 20;
Ny = 20;
Nz = 1;

G = cartGrid([Nx Ny Nz],[x_size y_size z_size]);
G = computeGeometry(G);

%% We need two rock structures for the dual-porosity model: one for the fracture
% system and one for the non-fractured rock matrix. Fractures usually will 
% have high permeability and low pore volume, while matrix usually will 
% have a higher porosity and decreased permeability. Therefore,
% interconnected fractures generally serves as a "highway" for fluid flow, 
% while the fluid transfer with the matrix will happen at a larger
% timescale
kf = 10000*milli*darcy;
phif = 0.01;

km = 1*milli*darcy;
phim = 0.1;

% Rock Structures
rock_fracture = makeRock(G, kf, phif);
rock_matrix = makeRock(G, km, phim);

%% Pressures of injector and producer wells
pres = 0*psia;
pprod = pres - 1000*psia;
pinj = pres + 1000*psia;

%% We also need two fluid structures for the dual-porosity model. Fractures are
% usually set as having zero capillary pressure, but we don't define any 
fluid_matrix = initSimpleADIFluid('mu', [1 2], 'rho', [1 1 1], 'n', [2 2 1]);
fluid_fracture = fluid_matrix;

bp = 1;
bi = 1.1;

m = (bi-bp)/(pinj-pprod);
b = bp-m*pprod;

fluid_fracture.bW = @(p)(b+m*p);
fluid_fracture.bO = @(p)(b+m*p);

%% Capillary Pressure
% Pe = 100*kilo*Pascal;
fluid_fracture.pcOW=@(swm)0;
fluid_matrix.pcOW=@(swm)-100*psia*swm+100*psia;

%% Here you can send information to the dual-porosity model that does not need to be calculated
% on each timestep
dp_info.constant_shape_factor = 50;

%% Add the wells
W = addWell([],G,rock_fracture,Nx*Ny,'type','bhp','Val',pprod,'Comp_i',[1,1]);
W = addWell(W,G,rock_fracture,1*1,'type','bhp','Val',pinj,'Comp_i',[1,0]);

%% Create the model. TwoPhaseOilWaterDPModel is a modified version of 
% TwoPhaseOilWaterModel that adds dual porosity behaviour
gravity off

%% Model definition
model = TwoPhaseOilWaterDPModel(G, rock_fracture, fluid_fracture,...
                               rock_matrix, fluid_matrix, []);
                               
% The shape factor and transfer function
fracture_spacing = repmat([1,1,1],G.cells.num,1);
shape_factor_name = 'KazemiShapeFactor';
model.transfer_model_object = GilmanKazemiOilWaterTransferFunction(shape_factor_name,fracture_spacing);

%% Initialize the field with a constant pressure and fully saturated by oil
state0.pressure = ones(G.cells.num,1)*pres;
state0.s = repmat([0 1],G.cells.num,1);
state0.swm = zeros(G.cells.num,1);
state0.pom = ones(G.cells.num,1)*pres;

%% Initialize the well solution
state0.wellSol= initWellSolAD(W, model, state0);
state = state0;

solver = NonLinearSolver();

%% Handles to pictures that are going to be plotted at each timestep
fig1 = figure('Position',[100,0,600,800]);
fig2 = figure('Position',[700,0,600,800]);

%% Source
src_val = 0.0001*sum(poreVolume(G,rock_fracture))/day;
src = addSource([], 200, src_val, 'sat', [0,1]);

%% Simulate the models
dT = 0.1*day;
n = 20;
for i = 1:n
    %% Advancing fields
    state = solver.solveTimestep(state, dT, model, 'W', W, 'src', src);
    disp(i)

    %% Plotting fields
    figure(fig1)
    subplot(2,1,1)
    title(['Saturation of Oil in the Fractures @ step = ' num2str(i)])
    p = plotCellData(G,1-state.s(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar
%     caxis([0 1])

    figure(fig1)
    subplot(2,1,2)
    title(['Saturation of Oil in the Matrix @ step = ' num2str(i)])
    p = plotCellData(G,1-state.swm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar
%     caxis([0 1])

    figure(fig2)
    subplot(2,1,1)
    title(['Transfer Rate of Water @ step = ' num2str(i)])
    p = plotCellData(G,state.Twm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar

    figure(fig2)
    subplot(2,1,2)
    title(['Transfer Rate of Oil @ step = ' num2str(i)])
    p = plotCellData(G,state.Tom(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)
    colorbar

    drawnow;
end
