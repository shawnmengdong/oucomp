function [problem, state] = pressureEquationOilWaterDP(state0, state, model, dt, drivingForces, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'propsPressure', [], ...
             'staticWells',  false, ...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

% assert(isempty(drivingForces.bc) && isempty(drivingForces.src))

s = model.operators;
f = model.fluid;

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellsol');

% Matrix properties
[pom,swm] = model.getProps(state, 'pom','swm');
[pom0,swm0] = model.getProps(state0, 'pom','swm');

%Initialization of independent variables ----------------------------------
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p,pom, wellVars{:}] = ...
            initVariablesADI(p, pom, wellVars{:});
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'pressure', 'pom', wellVarNames{:}};

p_prop = opt.propsPressure;
otherPropPressure = ~isempty(p_prop);
if ~otherPropPressure
    p_prop = p;
end

% -------------------------------------------------------------------------
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p_prop, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p_prop, sW, krW, T, gdz);
bW0 = f.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p_prop, sO, krO, T, gdz);
bO0 = getbO_BO(model, p0);

%% Properties for Matrix
pvMultm = pvMult;
pvMultm0 = pvMult0;

% Using capillary pressure information
pcOWm = 0;
pcOWm0 = 0;
if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
    pcOWm  = model.fluid_matrix.pcOW(swm);
    pcOWm0  = model.fluid_matrix.pcOW(swm0);
end

pwm = pom - pcOWm;
pwm0 = pom0 - pcOWm0;

% SMALL TO DO HERE: WE USE THE REL PERMS OF THE MATRIX TO EVALUATE THE
% EFFECTIVE PERMEABILITY
som = 1-swm;
som0 = 1-swm0;
[krWm, krOm] = model.evaluateRelPerm({swm, som});

bWm = f.bW(pwm);
bOm = f.bO(pom);

bWm0 = f.bW(pwm0);
bOm0 = f.bO(pom0);

%% Transfer
vb = model.G.cells.volumes;

matrix_fields.pom = pom;
matrix_fields.swm = swm;
fracture_fields.pof = p;
fracture_fields.swf = sW;

transfer_model = model.transfer_model_object;

[Talpha] = transfer_model.calculate_transfer(model,fracture_fields,matrix_fields);

Twm = vb.*Talpha{1};
Tom = vb.*Talpha{2};

%% Go on
if otherPropPressure
    % We have used a different pressure for property evaluation, undo the
    % effects of this on the fluxes.
    dp_diff = s.Grad(p) - s.Grad(p_prop);
    
    vW = -s.faceUpstr(upcw, mobW).*s.T.*(dpW + dp_diff);
    vO = -s.faceUpstr(upco, mobO).*s.T.*(dpO + dp_diff);
end

% These are needed in transport solver, so we output them regardless of
% any flags set in the model.
state = model.storeFluxes(state, vW, vO, []);
state = model.storeUpstreamIndices(state, upcw, upco, []);
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
	state.Twm = double(Twm);
    state.Tom = double(Tom);
end
% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% oil fracture:
oil_fracture = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + s.Div(bOvO);
oil_fracture = oil_fracture + Tom;

% water fracture:
wat_fracture = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);
wat_fracture = wat_fracture + Twm;

% oil matrix
oil_matrix = (s.pv_matrix/dt).*( pvMultm.*bOm.*som - pvMultm0.*bOm0.*som0 );
oil_matrix = oil_matrix - Tom;

% water matrix 
wat_matrix = (s.pv_matrix/dt).*( pvMultm.*bWm.*swm - pvMultm0.*bWm0.*swm0 );
wat_matrix = wat_matrix - Twm;

eqTmp = {wat_fracture, oil_fracture};
names = {'water', 'oil'};
types = {'cell', 'cell'};
rho = {rhoW, rhoO};
mob = {mobW, mobO};
sat = {sW, sO};

[eqTmp, state] = addBoundaryConditionsAndSources(model, eqTmp, names, types, state, ...
                                                                 {pW, p}, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

                                   
dissolved = {};
[eqTmp, names, types, state.wellSol] = model.insertWellEquations(eqTmp, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {}, dt, opt);

wat_fracture = eqTmp{1};
oil_fracture = eqTmp{2};

[eqs, names, types] = deal({});



eqs{1} = (dt./s.pv).*(oil_fracture./bO + wat_fracture./bW);
names{1} = 'pressure';
types{1} = 'cell';

eqs{2} = (dt./s.pv_matrix).*(oil_matrix./bOm + wat_matrix./bWm);
names{2} = 'pressure_matrix';
types{2} = 'cell';

state.timestep = dt;
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);


end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
