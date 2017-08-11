% test case
% liquid-vapor phase equilibia calculation
% Define the components and load pure physsical properties
z=[.01210 .01940 .65990 .08690 .05910 .09670 .04745 .01515 .00330];
% Define the thermodynamic models
T= 300; % [K]  temperature
p= 10e5; % [Pa] pressure
thermo = addThermo();
thermo.EOS = @PREOS; 
mixture = addMixture();
mixture.mole_fraction = z;
mixture.pressure = p;
mixture.temperature = T;
% Define flash options
options.convergence_eps = 1e-12;   %convergence tolerance for fugacity
options.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
options.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
options.max_outer_loop = 1000;   %max number of fugacity updates
% flash for surface properties
[success_flag,stability_flag,vapor_y_o,liquid_x_o,vapor_frac,zc,phase_index,cubic_time,b]=GI_flash(mixture,thermo,options);

