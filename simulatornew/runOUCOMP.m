%OUComp
%Load AD module
mrstModule add ad-fi

%Define Fluid Section
fluid = setupFluid();

%Define Grid Section
G = setupGrid();

%Define Rock Section
rock = setupRock(G);

%Other Parameters
param = setupParam();

%Define Wells
W = setupWells(G,rock);
