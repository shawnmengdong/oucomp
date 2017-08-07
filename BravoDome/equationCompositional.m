%%
%
% perlOneLineDescription(Residual equations functions)
%
% In this function, the residual equations are assembled and evaluated. By using
% _Automatic Differentiation_ (AD), we compute automatically the derivatives of the
% residual. Examples of use of automatic differentiation can be found in
% <http://www.sintef.no/Projectweb/MRST/Publications/ MRST book>.
%
% The governing equations are the mass conservation equations for the 4 components: He,
% Ne, CO2 and water, and the equilibrium equations, or flash equations.
%
% The mass conservation equations are discretized in space using a Two Point Flux
% Approximation. The transmissibilities are computed before hand when the system is set up
% in the function |setupSystem|. The equations are solved implicitly in time. 
%
%% 

function eqs = equationCompositional(state0, state, dt, bc, system, varargin)

    opt = struct('Verbose',     mrstVerbose,...
                 'scaling',     [],...
                 'history',     [],  ...
                 'iteration',   -1,  ...
                 'stepOptions', []);

    opt = merge_options(opt, varargin{:});

    s = system.s;
    f = system.fluid;
    Temp = system.Temp;
    T = s.T; % Transmissibilities.
    nComp = system.nComp;
    G = system.G;
    
    %% 
    % 
    % The main variables are pressure, liquid saturation and total concentrations. They
    % are converted to AD class objects by the function |initVariablesADI|.  Since they
    % are initialized as main variables, their Jacobian are identity matrices. By operator
    % overloading, any expression which contains these variables (|p|, |sL| and |C|)
    % will also be an AD object and have a derivative. In this way, the derivatives of
    % the residuals are automatically obtained.

    p    = state.pressure;
    sL   = state.sL;
    C    = state.C;

    [p, C{:}, sL] = initVariablesADI(p, C{:}, sL);

    p0   = state0.pressure;
    sL0  = state0.sL;
    C0   = state0.C;
    % The previous values should be stored in state
    [cg0, cl0] = computeComposition(C0, p0, sL0, system); 
    [Cw0, cwg0, cwl0] = computeWaterComp(p0, C0, cl0, sL0, system);
    

    %% 
    % Compute the relative permeabilities. We also need the values at the boundaries
    % where Dirichlet conditions hold.
    %
    
    [krL, krG] = f.relPerm(sL); % Liquid, gas relative permeability
    bd = bc.dirichlet; % short-cut for dirichlet bc values.
    [bc_krL, bc_krG] = f.relPerm(bd.s(:, 2)); % Liquid and gas relative permeability as
                                              % the boundary
    
    
    %%
    % Set up the gravity term. The molar mass are needed.
    %
    
    g  = 9.8;
    dz = s.dz;
    mmW   = f.mmW; % f.mmW returns molar mass for water.
    mmC   = f.mmC; % f.mmC returns a vector of molar mass for each component.

    
    %%
    % Compute the mobilities
    %
    
    mobL   = krL./f.muL;
    mobG   = krG./f.muG;
    bc_mobL   = bc_krL./f.muL;
    bc_mobG   = bc_krG./f.muG;

    %%
    % Compute the upstream direction for each component.
    %
    
    dpC = cell(nComp, 1);
    upC = cell(nComp, 1);
    for ic = 1:nComp
       dpC{ic} = s.p_grad(p) - g*(mmC(ic).*dz);
       upC{ic} = (double(dpC{ic})>=0);
    end
    dpW = s.p_grad(p) - g*(mmW.*dz);
    upW  = (double(dpW)>=0);
    
    
    %% 
    % Components the composition of the phases: the gas |cg| and liquid |cl|
    % concentrations. Note that |cg| and |cl| are AD objects and contain the Jacobian
    % with respect to the main AD variables |C|, |p| and |sL|.
    %
    
    [cg, cl] = computeComposition(C, p, sL, system);
    fluxC = cell(nComp, 1);
    
    %%
    % Compute the residual of the mass conservation equation for each component (He, Ne and
    % CO2)
    %
    
    for ic = 1 : nComp
       %%
       % The function |s.faceConcentrations| computes concentration on the faces given
       % cell concentrations, using upwind directions.
       %
       bc_val = bd.cg{ic}.*bc_mobG + bd.cl{ic}.*bc_mobL; 
       fluxC{ic} = s.faceConcentrations(upC{ic}, cg{ic}.*mobG + cl{ic}.*mobL, bc_val);
       eqs{ic} = (s.pv/dt).*(C{ic} - C0{ic}) + s.div(fluxC{ic}.*T.*dpC{ic});
    end
    
    %%
    % Compute the residual of the mass conservation equation for water.
    %
    
    [Cw, cwg, cwl] = computeWaterComp(p, C, cl, sL, system);
    bc_val = bd.cw(:, 1).*bc_mobG + bd.cw(:, 2).*bc_mobL;
    fluxW = s.faceConcentrations(upW, cwg.*mobG + cwl.*mobL, bc_val);
    eqs{nComp + 1} = (s.pv/dt).*(Cw - Cw0) + s.div(fluxW.*T.*dpW);

    %%
    % Compute the residual for the flash equation for the saturation. Note that, since the
    % function |equationCompositional| is always called after the flash are equations are
    % solved, this residual is always equal to zero in every cells. *But*, the residual
    % variable is an AD objects which also contains the derivative. For this particular
    % residual, we obtain in this way, given in an implicit form, the derivative of the
    % liquid saturation |sL| with respect to the total concentrations |C| and the
    % pressure |p|.
    %
    
    vp = system.vp;
    eqs{nComp + 2} = (1 - p/vp);
    R = system.R;
    k = system.k;
    for ic = 1:nComp
       eqs{nComp + 2} = eqs{nComp + 2} + 1/vp*R*Temp.*C{ic}./(1 + (R*Temp/k(ic) - 1).*sL);
    end
    
    %% 
    % Check the status for each cell
    %
    % * |st = 1| : Pure liquid phase
    % * |st = 2| : Pure gas phase
    % * |st = 3| : Both liquid and gas phase
    %
    
    st   = getCellStatus(state, system);
    
    %%
    % We do not handle the case only gas is present.
    %
    
    assert(all(st ~= 2), 'pure gas phase case is not handled');

    %%
    % For the cells where only liquid is present (|st == 1|), the liquid saturation
    % equation is simly |sL = 1|.
    %
    
    is_st_one = (st == 1);
    if any(is_st_one)
       eqs{nComp + 2}(is_st_one) = sL(is_st_one) - 1;
    end
    
    %%
    % Add input fluxes. boundary fluxes are handled here
    % 
    
    for ic = 1 : nComp
       eqs{ic}(bc.influx_cells) = eqs{ic}(bc.influx_cells) - bc.C_influx{ic};
    end
    eqs{nComp + 1}(bc.influx_cells) = eqs{nComp + 1}(bc.influx_cells) - bc.water_influx;
    
end

function st = getCellStatus(state, system)
   
%%
% This function computes the status in each cell
%
% * |st = 1| : Pure liquid phase
% * |st = 2| : Pure gas phase
% * |st = 3| : Both liquid and gas phase
%
   
   [R, k, Temp, vp, nComp] = deal(system.R, system.k, system.Temp, system.vp, system.nComp);
   p = state.pressure;
   
   %%
   % Initiate status by assuming that gas and liquid are both present.
   %
   
   st = 3*ones(numel(p), 1);

   %%
   % Set pure liquid status, if pressure is above threshold.
   % 
   
   C = cell2mat({state.C{:}});
   st(C*k + vp <= p) = 1;

   %%
   % Set pure liquid status, if pressure is below threshold.
   % 
   
   omega_g = R*Temp*ones(nComp, 1);
   min_pressure = C*omega_g + vp;
   st(p < min_pressure) = 2;
   
end

