% The main state variables are :
% * the total molar fraction for each component(|state.z|)
% * the water satruation (|state.sw|)
% * the pressure (|state.p|),
% * the total molar density(|state.F|)
% There are in total of Nc+3 variables and equations
%
%%
% Load automatic differentiation module.
run('OUCOMPInitialize.m'); %first add all paths
mrstModule add ad-fi 

%%
% Production parameters. 
param.outflux_p = 4e6; %production pressure

%Reservoir initial conditions
param.Temp = 30 + 273.15; % Temperature (in Kelvin)
param.z_initial = [50, 50]/100;  % initial Molar Fraction (CH4,C8H18)
param.sw_initial = 0.25;
param.so_initial = 0.5;
param.sg_initial = 0.25;
param.p_initial = 10e6;  %initial reservoir pressure

% Time-step parameters
param.dt = 200*day;             % Time step
param.total_time = 100000*day;  % Total time

% Parameters for the non-linear solver
param.maxIterations = 50;

%% 
% create variable |G| and |rock|
[nx,ny,nz] = deal( 10, 10, 10);
[Lx,Ly,Lz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz ]);
G = computeGeometry(G);
rock = makeRock(G, 30*milli*darcy, 0.3);


%%
% Set up fluids
components_formula = {'CH4','C10H22'}; %main component except for water
fluid = setupFluids(components_formula);

%%
% Set up wells
gravity = 9.8; %gravity term
W = setupWells(G,rock);

%%
% Set up systems 
system = setupSystem(G,rock,param,fluid,W);

%%
% Setup the initial state.
state0 = inistate(system);

total_time = param.total_time;
dt         = param.dt;
steps      = dt*ones(floor(total_time/dt), 1);
t          = cumsum(steps);

%% 
% Time step iterations.
for tstep = 1 : numel(steps)
   dt = steps(tstep);
   % Call non-linear solver solvefi
   [state, conv] = solvefi(state0, dt, system, @equationCompositional);
   if ~(conv)
      error('Convergence failed. Try smaller time steps.')
   end
   state0 = state;
end



 
