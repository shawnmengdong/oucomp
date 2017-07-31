function system = setupSystem(G,rock,param,fluid,W)
% Set system variables
system.G         = G; %grid variable
system.fluid     = fluid;
system.nComp     = fluid.num_mixComponent; % 3 components (CO2£¬CH4,C8H18) except for water
system.cellwise  = 1:system.nComp+3; % number of equations used
system.param = param;
system.nonlinear = setNonlinearSolverParameters(param);
system.podbasis  = [];
system.Temp      = param.Temp*ones(G.cells.num,1);
%Set static flash options
system.thermo = addThermo();
system.thermo.EOS = @PREOS; %default EOS is Peng-Robinson
system.mixture = addMixture(fluid.components);
system.flashopt.convergence_eps = 1e-12;   %convergence tolerance for fugacity
system.flashopt.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
system.flashopt.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
system.flashopt.max_outer_loop = 1000;   %max number of fugacity updates
system.W = W;
%Set system functions
system.s = sysFunctions(G, rock,W);

end