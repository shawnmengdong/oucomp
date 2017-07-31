% test case
% liquid-vapor phase equilibia calculation
% Define the components and load pure physsical properties
z=[0.2,0.4,0.4];
z=mynormalize(z);
[component, comp_flag] = addComponents({'CO2','CH4','C8H18','H2O'});
% Define the thermodynamic models
T= 300; % [K]  temperature
p= 10e5; % [Pa] pressure
thermo = addThermo();
thermo.EOS = @PREOS; 
mixture = addMixture(component, T, p);
mixture.mole_fraction = z;
mixutre.pressure = p;
mixture.temperature = T;
% Define flash options
options.convergence_eps = 1e-12;   %convergence tolerance for fugacity
options.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
options.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
options.max_outer_loop = 1000;   %max number of fugacity updates
% flash for surface properties
[success_flag,stability_flag,vapor_y,liquid_x,vapor_frac,cubic_time]=GI_flash(mixture,thermo,options);

