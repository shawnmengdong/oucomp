%Initialize the state variables

function state = inistate(system)
   param = system.param;
   G = system.G;
   W = system.W;
   nComp = system.nComp; % Number of components (without water): 3 (CO2,CH4,C8H18)
   ini_pressure = param.p_initial; %initial reservoir pressure
   nc = G.cells.num;
   Temp = param.Temp*ones(nc,1);
   pressure = ones(nc,1)*ini_pressure;
   z = cell(1,nComp);
   
   for ic = 1 : nComp
      z{ic} = param.z_initial(ic)*ones(nc, 1);                 
   end
   
   [~,~,~,molar_density,~,~] = flash_calculation(pressure,Temp,z,system);

   state.z = z;
   state.pressure = pressure;
   state.sw = ones(nc, 1)*param.sw_initial;
   state.F = param.so_initial*molar_density(:,1)+param.sg_initial*molar_density(:,2);
   state.bhp = W(1).val;
   state.q = 0;
   
end