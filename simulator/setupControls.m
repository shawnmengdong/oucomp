%%
% There are two types of controls:
%
% * Dirichlet boundary conditions, where the pressure is given on selected external faces.
% * Neumann boundary conditions, where sources/sinks are included at given cells. We have thus not
% included well models.
%
% The controls locations have been hard coded here, but they can be changed. We have
% chosen a set of four cells on the western side for the injection and a face on the
% eastern side where we impose pressure.
%%

function bc = setupControls(G,param,fluid)
   

   nComp = fluid.num_mixComponent; % 3 components (CO2,CH4,C8H18) except for H2O
   outflux_p = param.outflux_p; %outflux pressure
   influx_p = param.influx_p;
   influx_z = param.influx_z;  %injection molar fraction
   influx_rate =  param.influx_rate;
   Temp = param.Temp;
   %% outflux boundary condition
   % Given pressure on eastern face of the cell 411, which is the most eastern cell of
   % the reservoir.
   
   east_cell = 411;
   east_cell_faces = G.cells.faces(G.cells.facePos(east_cell):G.cells.facePos(east_cell + 1) -1 , :);
   bc.dirichlet.faces =  east_cell_faces(east_cell_faces(:, 2) == 2 , 1);
   bc.dirichlet.pressure = outflux_p;
   bc.dirichlet.cell =east_cell;

   %%
   % Given influx on a set of 4 cells on the western side of the reservoir.
   %
   west_cells = [4527; 4528; 4451; 4452];
   nwc = numel(west_cells);
   bc.influx_cells = west_cells;
   % Since in this case we are only injecting CO2, it is single gas phase
   % with yi = zi
   influx_z = repmat(influx_z', nwc, 1);
   cell_volume = G.cells.volumes(west_cells);   %in m^3
   
   for ic = 1:nComp
      bc.Component_influx{ic}  =  influx_p*influx_rate/8.314/Temp./cell_volume.*influx_z(:, ic); %in mol/m3/s
   end
   bc.Total_influx = influx_p*influx_rate/8.314/Temp./cell_volume; %in mol/m3/s
   
   %bc.water_influx = param.influx_rate*Cw; 
   
end



