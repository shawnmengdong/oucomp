% Setup fluid properties
%
function fluid = setupFluids(components_formula)










% We consider a quadratic relative permeability curve
fluid.relPerm = @(sL) quadraticRelPerm(sL);
% Defining fluid molar mass
num_mixComponent = length(components_formula);
componentC = addComponents(components_formula); %component structure for components except for water
componentW = addComponents({'H2O'});
fluid.mmW = componentW.MW;  % molar mass of H20
mmc = zeros(1,length(components_formula));
for i = 1:num_mixComponent
    mmc(i) = componentC(i).MW;
end
fluid.mmC = mmc; % molar masses of CO2,CH4 and C8H18 in kg/mol

% Fluid viscosity for the liquid and gas (Assuming constant for now)
fluid.muL = 1e-3;  % Liquid
fluid.muW = 1e-3;  % Water  Pa.S
fluid.muG = 1e-5;  % Gas

fluid.epw = 5.55e4; %molar density of water, 55000mol/m3
fluid.rhow = 1000; %1000kg/m3

litre = 1e-3*meter^3;
rho = 1*kilogram/litre;
fluid.mv = fluid.mmW/rho;
fluid.num_mixComponent = num_mixComponent;
fluid.components=componentC;
fluid.nPhase_max = 2; %except for water phase
end