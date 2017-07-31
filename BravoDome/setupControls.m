%%
%
% There are two types of controls:
%
% * Dirichlet boundary conditions, where the total concentrations and pressure is given
% on selected external faces.
% * Injection fluxes, where sources are included at given cells. We have thus not
% included well models.
%
% The controls locations have been hard coded here, but they can be changed. We have
% chosen a set of four cells on the western side for the injection and a face on the
% eastern side where we impose pressure.
%%

function bc = setupControls(G, fluid, param)
   
   system.R         = getR();
   system.k         = getHenryCoef();
   system.Temp      = param.Temp;
   system.vp        = vaporPressure(system.Temp); 
   system.fluid     = fluid;
   nComp = 3; % 3 components (He, Ne, CO2)
   system.nComp     = nComp;
   
   influx_p = param.influx_p;
   outflux_p = param.outflux_p;
   
   %%
   % We assume injection fluid is only composed of (ideal) gas
   %
   
   influx_C = param.influx_C*influx_p/(system.R*system.Temp);

   %% 
   % We assume the output fluid is only composed of liquid (typically same composition as
   % initial state). This assumption should actually not be used if there is only
   % outflow on this faces. Indeed, in this case, the boundary values of the
   % concentrations are not used due to upwinding face evaluation. 
   %
   
   omega = cell2mat(omega_l(outflux_p, system));
   outflux_C = 1/(omega'*param.outflux_C)*param.outflux_C;
   outflux_C = outflux_C(1:nComp);

   %%
   % Given pressure on eastern face of the cell 411, which is the most eastern cell of
   % the reservoir.
   %
   
   east_cell = 411;
   east_faces = G.cells.faces(G.cells.facePos(east_cell):G.cells.facePos(east_cell + 1) -1 , :);
   bc.dirichlet.faces =  east_faces(east_faces(:, 2) == 2 , 1);
   bc.dirichlet.pressure = outflux_p;
   bc.dirichlet.C = arrayfun(@(x)(x), outflux_C, 'uniformoutput', false); 
   [cg, cl, cw, s, Cw] = flash_calculation(outflux_C', outflux_p, system);
   bc.dirichlet.cg = arrayfun(@(x)(x), cg, 'uniformoutput', false); 
   bc.dirichlet.cl = arrayfun(@(x)(x), cl, 'uniformoutput', false); 
   bc.dirichlet.cw = cw; 
   bc.dirichlet.s  = s; 
   bc.dirichlet.Cw = Cw; 

   %%
   % Given influx on a set of 4 cells on the western side of the reservoir.
   %
   
   west_cells = [4527; 4528; 4451; 4452];
   nwc = numel(west_cells);
   bc.influx_cells = west_cells;
   influx_C = repmat(influx_C', nwc, 1);
   influx_p = repmat(influx_p, nwc, 1);
   [cg, cl, cw, s, Cw] = flash_calculation(influx_C, influx_p, system);
   bc.dirichlet.C  = cell(nComp, 1);
   for ic = 1:nComp
      bc.C_influx{ic}  = param.influx_rate*influx_C(:, ic);    %There is a problem with it, fix it£¡
   end
   bc.water_influx = param.influx_rate*Cw; 
   
end



