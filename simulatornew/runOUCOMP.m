%OUComp
%Load AD module
mrstModule add ad-fi

%Define Fluid Section
fluid = setupFluid();

%Define Grid Section
G = setupGrid();

%Define Rock Section
rock = setupRock(G);

%Other Parameters
param = setupParam();

%Define Wells
W = setupWells(G,rock);

%Define system
system = setupSystem(G,rock,param,fluid,W);

%Initialize first state
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
