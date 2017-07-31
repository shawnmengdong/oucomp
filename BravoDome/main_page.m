%% 
% 
% 
% The Bravo Dome carbon dioxide gas field is located in northeast New Mexico. The Bravo Dome field
% covers approximately 3200 km^2. Production in 1989 was 3.2e9 cubic meter of gas from 272
% wells. Cumulative production at the end of 1989 was 17e9 cubic meter. Estimated recoverable
% reserves are more than 283e9 cubic meter. The gas is 98-99% CO2. Most CO2 produced from Bravo Dome
% is used for enhanced oil recovery in the Permian basin. We thank
% <http://www.jsg.utexas.edu/student/kiran_sathaye/ Kiran Sathaye> for his help in setting up this
% tutorial.
%
%
% We consider a compositional system composed of Helium, Neon, Carbon dioxide and
% Water. We use Henry's law to compute the ratio of each component in the gas and liquid
% phase. For water, we assume that the water vapor pressure is constant. We refer to the
% note <perlAddNoteLink() Bravo tutorial note> for a description of the model with the
% assumptions that are used.
%
% The grid is constructed from raw data tables consisting of arrays (dimension 100x100)
% which contains the permeability, the porosity, the depths of the lower and upper layers
% for a Cartesian grid. See function perlAddLink(setupGeometry).
%
% <html><table><tr><td><img src="permeability.png" width="500"/></td><td><img
% src="porosity.png" width="500" /></td></tr></table></html>
%
% Below, we list the files that are needed to set up this compositional solver. In the
% first list, we have the main functions that have been written for this tutorial. In the
% list _helper functions_, we find short functions that could have been included in the
% _main functions_ but which , for the sake of clarity, have been written as separate
% functions. At last, we list the functions from MRST (core and ad-fi modules) that are
% used in this tutorial. We use a fully-implicit in time and two point flux in space
% discretization.
%
% We observe that
%
% * The size of the code needed to implement this compositional simulator on top of
% MRST is particularly short.
%
% * Moreover, only a few functions from MRST are in fact used. They reflect two of the
% most interesting features of MRST when it comes to rapid prototyping: a flexible and
% robust grid structure and a set of handy AD tools.
% 
%
% <html><h2>Main functions</h2></html>
%
% * perlAddLink(runBravo) : In this function, the parameters for the computation are
% set. The geometry is loaded. The fluid properties are set. The non-dynamical
% computational structures (transmissibilities, discrete differential operators ) are
% computed, by calling perlAddLink(setupSystem). The control variables (injection rate and
% composition, output pressure) are set up by calling perlAddLink(setupControls). The
% schedule is set up. Finally, the computation is started by solving at each time step the
% fully-implicit system.
%
% * perlAddLink(solvefi): This function is a one-step Newton solver. The equations are
% assembled using perlAddLink(equationCompositional).
%
% * perlAddLink(flash_calculation): This function computes the liquid saturation and the
% concentrations of all the components in the gas and liquid phase, given the total
% concentrations and the pressure.
%
% * perlAddLink(equationCompositional): The is the function where the equations are
% assembled. The equations are four mass conservation equations, one for each element (He,
% Ne, CO2, H2O), and one equilibrium equation, the saturation equation. The residual
% are evaluated. Using the _automatic differentiation_ framework, the Jacobians of the
% residual equations is also automatically computed. 
%
% * perlAddLink(setupControls):  In this function, we set the control variables: We
% define the cells where the injection is performed and the faces where the output
% pressure is set. The injection rate and output pressure value are set.
%
% * perlAddLink(setupGeometry): Set up the geometry from raw data.
%
% * perlAddLink(setupSystem): In this function, we compute the transmissibilities and set
% up the discrete differential operators |div| and |grad|.
%
% <html><h2>Helper functions</h2></html>
%
% * perlAddLink(computeComposition): Computes the concentrations in the gas and liquid
% phase of all elements except water, given the total concentration, pressure *and* liquid
% saturation
%
% * perlAddLink(computeWaterComp): Compute the concentrations of water in both phases
% given given the total concentration, pressure and liquid saturation.
%
% * perlAddLink(getHenryCoef): Returns the Henry solubility coefficients for each component.
%
% * perlAddLink(getR): Returns the ideal gas constant.
%
% * perlAddLink(initStateBravo): Initialize the state variables.
%
% * perlAddLink(omega_l): Computes the molar volume for the liquid phase.
%
% * perlAddLink(vaporPressure): Return the water vapor pressure.
%
% * perlAddLink(quadraticRelPerm): Set up quadratic relative permeabilities,
%
% * perlAddLink(setNonlinearSolverParameters): Set up nonlinear parameters for the solver.
%
%
% perlIncludeTopMenu()
%
% <html><h2>MRST functions used in this tutorial</h2></html>
%
% perlIncludeCoreFunctions() 
% 
% <html><h2>Download tutorial files</h2></html>
%
% <perlAddLinkToSource() tutorial-Bravo.tar.gz>
%
% <html><h2>Numerical results</h2></html>
% 
% 
% Initially, we set the initial molar fractions of neon and water equal to 0.5% and 99.5%
% respectively, with no helium and carbon dioxyde. The pressure in the whole reservoir is
% set to the output pressure, that is, 80 bar. We inject a mixture composed of 10%
% Helium and 90% CO2 at a rate of 1000m^3/day. To compute the injected mole amount, we
% solve the flash equations assuming an injection pressure of 130 bar.
%
% The injection is done on four cells on the western side of the reservoir and the output
% pressure is set on one face at the most eastern cell, see perlAddLink(setupControls).
%
% The simulation is run for 10000 days with a time step of 200 days. The following plots
% show the result of the simulation at end time. We observe that, as espected, the helium
% travels faster than the carbon dioxyde because of its lower solubility in the liquid
% phase.
%
% <html><table><tr><td><img src="plotC_He.png" width="500"/></td><td><img
% src="plotC_Ne.png" width="500"/></td><tr><td><img src="plotC_CO2.png" width="500"
% /></td><td><img src="plotC_sL.png" width="500"/></td></tr></table></html>
% 
%
%%