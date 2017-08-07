%%
%
% perlOneLineDescription(Single step non-linear solver)
% 
% Compute the next |state| of the system for given previous state |state0| and time step
% |dt|.
%
% The discretized residual equation are given by |equation| which is assembled in the
% function |equationCompositional|. 
%%

function [state, convergence] = solvefi(state0, dt, bc, system, equation, param,  varargin)

   opt = struct('verbose', false);
   opt = merge_options(opt, varargin{:});

   fluid = system.fluid;
   
   meta = struct('converged'  , false                       , ...
                 'stopped'    , false                       , ...
                 'relax'      , system.nonlinear.relaxation , ...
                 'stagnate'   , false                       , ...
                 'iteration'  , 0                           , ... 
                 'res_history', []                          );

   timer = tic;
   
   converged = false;
   state = state0;

   fprintf('%13s%-26s%-36s\n', '', 'CNV (oil, water)', 'MB (oil, water)');
   
   equation = @(state) equation(state0, state, dt, bc, system);
   
   %%
   % We start with the Newton iterations
   
   while ~ (converged || meta.stopped),
      % Save iteration number in meta info
      meta.iteration = meta.iteration + 1;

      %% 
      % Initial saturation solve
      % 
      % At each Newton step, we start by solving the flash equations and update the liquid
      % saturation variable.
      % 
      
      [C, p] = deal(state.C, state.pressure);
      init_sL = state.sL;
      Cvec = cell2mat(C);
      [~, ~, ~, s, ~] = flash_calculation(Cvec, p, system, init_sL);
      state.sL = s(:, 2);
      
      %%
      % The residual equations for the whole system (pressure, total concentrations,
      % liquid saturation) are assembled.
      %
      
      eqs = equation(state);

      %%
      % We call a standard linear solve to compute the Newton step |dx|.
      
      dx = SolveEqsADI(eqs, []);
      
      %%
      % We update |state|.
      %
      state      = updateState(state, dx, system);
      
      %%
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

%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%

function state = updateState(state, dx, system)

   nComp = system.nComp;
   dp = dx{1};
   for ic = 1 : nComp
      dC{ic} = dx{ic + 1};
   end
   dsL = dx{nComp + 2};
    
   step = 1;

   state.pressure = state.pressure + step*dp;
   for ic = 1 : nComp
      state.C{ic} = max(0, state.C{ic} + step*dC{ic});
   end
   state.sL = min(1, max(0, state.sL + step*dsL));
   
end
