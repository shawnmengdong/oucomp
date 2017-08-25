%% 
% Update function
%
% Given a Newton step, the state variables are updated. We cap the saturation variable.
%

function state = updateState(state, dx, system)

   nComp = system.nComp;
   dz = cell(nComp,1);
   for ic = 1 : nComp
      dz{ic} = dx{ic};
   end
   dp = dx{ic+1};
   dF = dx{ic+2};
   dsw = dx{ic+3};
   dbhp = dx{ic+4};
   dq = dx{ic+5};
   
   
   alfa = 1; %learning rate/update step size

   state.pressure = state.pressure + alfa*dp;
   state.F = state.F +alfa*dF;
   state.bhp = state.bhp + alfa*dbhp;
   state.q = state.q +alfa*dq;
   for ic = 1 : nComp
      state.z{ic} = min(1,max(0, state.z{ic} + alfa*dz{ic})); %0<zi<1
   end
   
   state.sw = min(1, max(0, state.sw + alfa*dsw));
   
end
