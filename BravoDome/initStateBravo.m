%%
%
% perlOneLineDescription(Initialize the state variables)
% 
%%


function state = initStateBravo(G, system, param)

   nComp = 3; % Number of components (without water): 3 (He, Ne, CO2)
   pressure = param.outflux_p;

   omega = cell2mat(omega_l(pressure, system));
   C_initial = 1/(param.C_initial*omega)*param.C_initial;
   C_initial = C_initial(1:nComp);
   
   [cg, cl, cw, s, Cw] = flash_calculation(C_initial, pressure, system);

   C = cell(nComp, 1);
   nc = G.cells.num;
   for ic = 1 : nComp
      state.C{ic} = C_initial(ic)*ones(nc, 1);                 
   end
   state.pressure = pressure*ones(nc, 1);
   state.sL = ones(nc, 1)*s(2);

end