%%
%
% Single step non-linear solver
% 
% Compute the next |state| of the system for given previous state |state0| and time step
% |dt|.
%
% The discretized residual equation are given by |equation| which is assembled in the
% function |equationCompositional|. 
%%

function [state, convergence] = solvefi(state0, dt, system, equation)

   param = system.param;
   meta = struct('converged'  , false                       , ...
                 'stopped'    , false                       , ...
                 'relax'      , system.nonlinear.relaxation , ...
                 'stagnate'   , false                       , ...
                 'iteration'  , 0                           , ... 
                 'res_history', []                          );

   timer = tic;
   converged = false;
   state = state0;
   equation = @(state) equation(state0, state, dt, system);
   
   % We start with the Newton iterations
   while ~ (converged || meta.stopped)
      % Save iteration number in meta info
      meta.iteration = meta.iteration + 1;
      % At each Newton step, we start by solving the flash equations 
      [pressure,Temp,comp] = get_flash_input(state,param);
      [liquid_x,vapor_y,L,molar_density,rho,~] = flash_calculation(pressure,Temp,comp,system);
      
      %update the PVT variable
      state.X = {liquid_x,vapor_y};
      state.phase_frac = [L,1-L];
      state.ep = molar_density;
      state.rho = rho;
      
      % The residual equations for the whole system
      eqs = equation(state);
      
      % We call a standard linear solve to compute the Newton step |dx|.
      dx = SolveEqsADI(eqs, []);
      
      % We update |state|.
      state      = updateState(state, dx, system);
      
      % We compute the residual values by calling |getResiduals|.
      % This function detects oscillation and stagnation
      
      [meta, residuals] = getResiduals(meta, eqs, system, false);

      %%
      % We test for convergence. Here we use a very stringent test given by a max norm
      %
      fprintf('Newton iteration %d, max residual: %5g\n', meta.iteration, ...
              max(abs(residuals)));
      converged = all(residuals <= system.nonlinear.tol); 

      if meta.stagnate,
         warning('newt:stagnate', 'Non-linear solver stagnated...')
      end
      
      if meta.iteration > system.nonlinear.maxIterations
         meta.stopped = true;
      end
      
   end


   if meta.stopped,
      warning('newt:maxit', ...
              ['Non-linear solver did not converge, stopped ', ...
               'by max iterations...']);
      convergence = false;
   else
      convergence = true;
   end
   
   dispif(mrstVerbose, 'Completed %d iterations in %1.2f s\n', ...
          meta.iteration, toc(timer));

end


