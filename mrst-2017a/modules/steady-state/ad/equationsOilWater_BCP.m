function [problem, state] = equationsOilWater_BCP(state0, state, model, dt, drivingForces, varargin)
% Get linearized problem for oil/water system with black oil-style
% properties
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
bcp = drivingForces.bcp;
s = model.operators;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0] = model.getProps(state0, 'pressure', 'water');

pBH    = vertcat(wellSol.bhp);
qWs    = vertcat(wellSol.qWs);
qOs    = vertcat(wellSol.qOs);

% Initialize independent variables.
if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
    else
        zw = zeros(size(pBH));
        [p0, sW0, zw, zw, zw] = initVariablesADI(p0, sW0, zw, zw, zw); %#ok
        clear zw
    end
end
% We will solve for pressure, water saturation (oil saturation follows via
% the definition of saturations) and well rates + bhp.
primaryVars = {'pressure', 'sW', 'qWs', 'qOs', 'bhp'};

% Evaluate relative permeability
sO  = 1 - sW;
sO0 = 1 - sW0;

[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity contribution
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO_BCP(model, p, p, sW, krW, T, gdz, bcp);
bW0 = model.fluid.bW(p0);

% Evaluate oil properties
[vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO_BCP(model, p, p, sO, krO, T, gdz, [], [], bcp);
bO0 = getbO_BO(model, p0);

if model.outputFluxes
    state = model.storeFluxes(state, vW, vO, []);
end
if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
    state = model.storeUpstreamIndices(state, upcw, upco, []);
end

% EQUATIONS ---------------------------------------------------------------
% Upstream weight b factors and multiply by interface fluxes to obtain the
% fluxes at standard conditions.
bOvO = s.faceUpstr(upco, bO).*vO;
bWvW = s.faceUpstr(upcw, bW).*vW;

% Conservation of mass for water
water = (s.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + s.Div(bWvW);

% Conservation of mass for oil
oil = (s.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + s.Div(bOvO);

eqs = {water, oil};
names = {'water', 'oil'};
types = {'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
eqs = addFluxesFromSourcesAndBC(model, eqs, ...
                                       {pW, p},...
                                       {rhoW,     rhoO},...
                                       {mobW,     mobO}, ...
                                       {sW, sO}, ...
                                       drivingForces);
% Finally, add in and setup well equations
if ~isempty(W)
    wm = WellModel(W);
    if ~opt.reverseMode
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [model.fluid.rhoWS, model.fluid.rhoOS];
        bw   = {bW(wc), bO(wc)};
        mw   = {mobW(wc), mobO(wc)};
        s = {sW(wc), sO(wc)};

        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, W, wellSol, ...
                                             pBH, {qWs, qOs}, pw, rhos, bw, mw, s, {},...
                                             'nonlinearIteration', opt.iteration);
        % Store the well equations (relate well bottom hole pressures to
        % influx).
        eqs(3:4) = weqs;
        % Store the control equations (trivial equations ensuring that each
        % well will have values corresponding to the prescribed value)
        eqs{5} = ctrleqs;
        % Add source terms to the equations. Negative sign may be
        % surprising if one is used to source terms on the right hand side,
        % but this is the equations on residual form.
        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
        eqs{2}(wc) = eqs{2}(wc) - cqs{2};
        
        names(3:5) = {'waterWells', 'oilWells', 'closureWells'};
        types(3:5) = {'perf', 'perf', 'well'};
    else
        [eqs(3:5), names(3:5), types(3:5)] = wm.createReverseModeWellEquations(model, state0.wellSol, p0);
    end
end


% PERIODIC BC & INCOMPRESSIBLE FLUID --------------------------------------

% In the case of periodic boundary conditions AND incompressible fluid, we
% do as in incompTPFA and set a level for the pressure to avoid singular
% Jacobian matrix.
if isempty(W) && ~isempty(bcp)
    % We are dependent on user input. The field 'isIncomp' must be set to
    % true if the fluid is incompressible.
    if isfield(model.fluid, 'isIncomp') && model.fluid.isIncomp
        if eqs{1}.jac{1}(1) > 0
            eqs{1}.jac{1}(1) = 2*eqs{1}.jac{1}(1);
        else
            [~, j] = max(diag(eqs{1}.jac{1}));
            eqs{1}.jac{1}(j,j) = 2*eqs{1}.jac{1}(j,j);
        end
    end
end



problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
