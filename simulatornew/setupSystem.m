function system = setupSystem(G,rock,param,fluid,W)

system.G = G; %Grid Variable
system.fluid = fluid;
system.param = param;
system.W = W;
system.rock = rock;
%-------------Nonlinear option------------------
system.nonlinear.maxIterations = 50;
% Relaxation parameters for Newton iterations
system.nonlinear.relaxation  = 1;
system.nonlinear.relaxMax    = 0.2;
system.nonlinear.relaxInc    = 0.1;
system.nonlinear.relaxRelTol = 0.01; %.1;
% Parameters for Newton's iterations
system.nonlinear.linesearch = false;
system.nonlinear.relaxType = 'sor';
system.nonlinear.tol = 1.0e-10;
%-----------------------------------------------------
system.s = setupFunctions(system);

end