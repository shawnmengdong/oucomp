function obj = pressurePenalizer(model, states, schedule, penalty, plim, varargin)
% states.pressure is a cell array.
% schedule is only used for time steps.
% penalty is a scalar.
% plim is a cell array.


% format of objective function:
%   obj = max(0, sign(p - plim)) * penalty * (p - plim)^2
% obj is computed for each time step (for both inj and mig periods).


   opt.ComputePartials = false;
   opt.tStep = [];
   opt = merge_options(opt, varargin{:});
   
   num_timesteps = numel(schedule.step.val);
   tSteps = opt.tStep;
   if isempty(tSteps)
      numSteps = numel(states);
      tSteps = (1:numSteps)';
      dts = schedule.step.val;
   else
      assert(numel(tSteps) == 1);
      numSteps = 1;
      dts = schedule.step.val(opt.tStep);
   end
   
   obj = repmat({[]}, numSteps, 1);
   k   = 2;
   max_amount_surp      = zeros(numSteps, 1);
   max_amount_surp_cinx = zeros(numSteps, 1);
   min_amount_under     = zeros(numSteps, 1);
   min_amount_under_cinx= zeros(numSteps, 1);
   for step = 1:numSteps     
      state = states{tSteps(step)}; %@@ +1?      
      p = state.pressure;
      % keep track of amount over or amount under plim at each time step
      max_amount_surp(step) = max(0,max(p-plim));
      if max_amount_surp(step) > 0
          [~,cinx] = max(p-plim);
          max_amount_surp_cinx(step) = cinx;
      end
      min_amount_under(step) = max(0,min(plim-p));
      if min_amount_under(step) > 0
         [~,cinx] = min(max(0,(plim-p)));
         min_amount_under_cinx(step) = cinx;
      end
      if opt.ComputePartials
        sG = state.s(:,2);   % place holders
        sGmax = state.sGmax; % place holders
        nW = numel(schedule.control(1).W);
        pBHP = zeros(nW, 1); % place holders
        qGs = pBHP;          % place holders
        qWs = pBHP;          % place holders
        [p, ~, ~, ~, ~, ~] = initVariablesADI(p, sG, sGmax, qWs, qGs, pBHP); 
      end
      dt = dts(step);
      tmp = max(0, sign(p - plim)) .* penalty .* (p - plim).^k/1e12; % scaling to MPa
      tmp = tmp .* model.G.cells.volumes;
      obj{step} = sum( tmp )./sum(model.G.cells.volumes);
      obj{step} = obj{step} * dt;
      if (tSteps(step) == num_timesteps)
      % no need to compute another portion of the obj fun here
         if ~opt.ComputePartials
             msg1 = 0;
             msg2 = 0;
             if any(max_amount_surp > 0)
                 [val,tinx] = max(max_amount_surp);
                 cinx = max_amount_surp_cinx(tinx);
                 msg1 = 100*val/plim(cinx);
             else
                 [val,tinx] = min(min_amount_under);
                 cinx = min_amount_under_cinx(tinx);
                 msg2 = 100*val/plim(cinx);
             end
             fprintf('Surpassed Plimit by %f (percent) of Plimit.\n', msg1)
             fprintf('Approached Plimit by %f (percent) of Plimit.\n', msg2)
         end
      end

   end
end