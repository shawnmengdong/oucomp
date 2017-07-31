%%
% 
% perlOneLineDescription()
%
% The main state variables are
%
% * the pressure (|state.p|),
% * the total concentrations of each component He, Ne and CO2 (|state.C|)
% 
% We assume that thermodynamical equilibrium is reached instantaneously. Hence, from the
% main state variables, we can obtain the other state variables (saturation,
% concentrations in gas and liquid phases) by solving the flash equations, see
% <perlAddNoteLink() Bravo tutorial note>.
%
%
% The flash equation can be reduced to a scalar equation in the liquid saturation
% variable. Thus, for the implementation convenience, we add 
%
% * the liquid saturation (|state.s|)
%
% to the main variable.
%%

%%
% Load automatic differentiation module.

mrstModule add ad-fi 


%% 
% Simulation parameters
% We store simulation parameters in a structure |param|

%%
% Injection parameter

param.influx_C  = [10; 0; 90]/100; 
param.influx_p  = 13e6; % for flash calculation, to compute influx water if inoutflux_unit='percent'.
param.influx_rate = 1000/day; % in m^3/s


%%
% Production parameters. 
% In this case only the pressure is in fact needed

param.outflux_C = [0; 0.25; 0; 99.75]/100; % Should not be used if upwinding occurs for all
                                           % components.
param.outflux_p = 8e6;

param.Temp = 30 + 273.15; % Temperature (in Kelvin)

param.C_initial = [0, 0.5, 0, 99.5]/100;  % initial concentrations (He, Ne, CO2, H20)
                                          % in percent

%%
% Time-step parameters
%
param.dt = 200*day;             % Time step
param.total_time = 100000*day;  % Total time


%%
% Parameters for the non-linear solver
%

param.maxIterations = 50;


%%
% Parameters to save data.
%

param.do_save = false;
param.output_dir = 'output'; % directory where state variables are saved.

%% 
% Load grid and rock parameters
% create variable |G| and |rock|

load('GEOMETRY');
load('ROCK');

%% 
% Setup fluid properties
% 
% We consider a quadratic relative permeability curve.

fluid.relPerm = @(sL) quadraticRelPerm(sL);

%%
% The following molar mass values are taken from
% <http://fr.wikipedia.org/wiki/Tableau_p%C3%A9riodique#mediaviewer/Fichier:Tableau_de_classification_p%C3%A9riodique_des_%C3%A9l%C3%A9ments.png
% Wikipedia>

mmH  = 1.00794*gram;  % molar mass of Hydrogen
mmO  = 15.9994*gram;  % molar mass of Oxygen
mmC  = 12.0107*gram;  % molar mass of Carbon
mmHe = 4.0026*gram;   % molar mass of Helium
mmNe = 20.1797*gram;  % molar mass of Neon

fluid.mmW = 2*mmH + mmO;  % molar mass of H20
fluid.mmC = [mmHe; mmNe; mmC + 2*mmO]; % molar masses of He, Ne amd CO2
   

%%
% Fluid viscosity for the liquid and gas

fluid.muL = 1e-3;  % Liquid
fluid.muG = 1e-5;  % Gas


%%
% Compressibility coefficient for the liquid phase with a reference pressure. We assume
% that the compressibility is independent on the composition

fluid.cl    = 4.4e-5/atm; % Compressibility
fluid.p_ref = 1*atm;      % Reference pressure

%% 
% Compute molar volume at standard condition (for pure water)

litre = 1e-3*meter^3;
rho = 1*kilogram/litre;
fluid.mv = fluid.mmW/rho;

%%
% We can include gravity.

gravity = false;

bc = setupControls(G, fluid, param);
 
%%
% Set system variables
   
system.G         = G;
system.fluid     = fluid;
system.nComp     = 3; % 3 components (He, Ne, CO2)
system.s         = setupSystem(G, rock, bc, param);
system.cellwise  = 1:5; % Used in function getResiduals which checks convergence.
   
   
%% 
% We have collected the parameter for the nonlinear solvers in the function
% |setNonlinearSolverParameters|

system.nonlinear = setNonlinearSolverParameters(param);
system.podbasis  = [];

system.R         = getR();
system.k         = getHenryCoef();
system.Temp      = param.Temp;
system.vp        = vaporPressure(system.Temp); 
   
%%
% Setup the initial state.
   
state0 = initStateBravo(G, system, param);

total_time = param.total_time;
dt         = param.dt;
steps      = dt*ones(floor(total_time/dt), 1);
t          = cumsum(steps);


%% 
% Time step iterations.

for tstep = 1 : numel(steps)

   dt = steps(tstep);
   
   %%
   % Call non-linear solver perlAddLink(solvefi)
   
   [state, conv] = solvefi(state0, dt, bc, system, @equationCompositional, param);

   if ~(conv)
      error('Convergence failed. Try smaller time steps.')
      return
   end

   if param.do_save
      save(fullfile(param.output_dir, sprintf('state%05d.mat', tstep)), 'state');
   end
   state0 = state;
   
end



 
