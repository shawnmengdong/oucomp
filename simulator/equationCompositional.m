% Residual equations functions
%
% In this function, the residual equations are assembled and evaluated. By using
% _Automatic Differentiation_ (AD), we compute automatically the derivatives of the
% residual. 
%
% The governing equations are the 9.26,9.27,9.28 and 9.29 in computational
% methods textbook.
%
% The mass conservation equations are discretized in space using a Two Point Flux
% Approximation. The transmissibilities are computed before hand when the system is set up
% in the function |setupSystem|. The equations are solved implicitly in time. 
%
function eqs = equationCompositional(state0, state, dt,system)

    s = system.s; % System functions
    fluid = system.fluid;
    T = s.T; % Transmissibilities.
    nComp = system.nComp; 
    G = system.G;
    W = system.W;
    wc = W(1).cells;
    nPhase_max = fluid.nPhase_max; %number of phases except for water
    % 
    % The main variables are converted to AD class objects by the function |initVariablesADI|.  Since they
    % are initialized as main variables, their Jacobian are identity matrices. By operator
    % overloading, any expression which contains these variables will also be an AD object and have a derivative. In this way, the derivatives of
    % the residuals are automatically obtained.
    p    = state.pressure;
    z   = state.z;
    F    = state.F;
    sw = state.sw;
    rho = state.rho;
    ep = state.ep;
    x = state.x;
    y = state.y;
    L = state.L;
    bhp = state.bhp;
    q = state.q;
    X = {x,y};
    phase_frac = [L,1-L];
    phase_flag = state.phase_flag;
    
    [z{:},p, F,sw,bhp,q] = initVariablesADI(z{:},p,F,sw,bhp,q);
    F0   = state0.F;
    z0  = state0.z;
    sw0 = state0.sw;
    
    % Compute the relative permeabilities. 
    [krL, krG, krW] = fluid.relPerm(sw); % Liquid, gas relative permeability

    % Set up the gravity term. The molar mass are needed.
    g  = norm(gravity);
    dz = s.grad(G.cells.centroids(:,3));
    rhow   = fluid.rhow; % returns mass density for water.
    epw = fluid.epw; %water molar density
    % Compute the mobilities
    mobComp   = {krL./fluid.muL,krG./fluid.muG};
    mobW   = krW./fluid.muW;
    
    % Compute the upstream direction of effective potential for each phase.
    dp_phase = cell(nPhase_max, 1);
    up_phase = cell(nPhase_max, 1);
    
    for ip = 1:nPhase_max
       dp_phase{ip} = s.grad(p) - g*(s.avg(rho(:,ip)).*dz);
       up_phase{ip} = (double(dp_phase{ip})>=0);
    end
    
    dpW = s.grad(p) - g*(rhow*dz);
    upW  = (double(dpW)>=0);
    
    %%
    % Compute the residual of the molar conservation equation for each component 
    for ic = 1 : nComp
        flux = zeros(size(T));
        for ip = 1:nPhase_max
            if phase_flag(ip)==1
               % The function |s.faceConcentrations| computes concentration on the faces given cell concentrations, using upwind directions.
               flux = flux + s.faceConcentrations(up_phase{ip},mobComp{ip}.*ep(:,ip).*X{ip}(:,ic)).*T.*dp_phase{ip}; %liquid flux:mobL*epL*x
            end
        end
         eqs{ic} = (s.pv/dt).*(F.*z{ic} - F0.*z0{ic}) - s.div(flux);
    end
    %%
    % Compute the residual of mass conservation equation for total component
     flux = zeros(size(T));
      for ip = 1:nPhase_max
            if phase_flag(ip)==1
               % The function |s.faceConcentrations| computes concentration on the faces given cell concentrations, using upwind directions.
               flux = flux + s.faceConcentrations(up_phase{ip},mobComp{ip}.*ep(:,ip)).*T.*dp_phase{ip}; %liquid flux:mobL*epL*x
            end
      end
       eqs{nComp+1} = (s.pv/dt).*(F - F0) - s.div(flux);
    
    
    %%
    % Compute the residual of the mass conservation equation for water.
    flux = s.faceConcentrations(upW, mobW.*epw).*T.*dpW;
    eqs{nComp + 2} = (s.pv*epw/dt).*(sw - sw0) - s.div(flux);

    %%
    % Compute the residual for saturation balance
      for ip = 1:nPhase_max
          S_total = zeros(size(p.val));
            if phase_flag(ip)==1
                S_total = S_total+ F.*phase_frac(:,ip)./ep(:,ip);
            end
      end
    eqs{nComp + 3} = S_total +sw-1;

    %%
    % Compute the residual for well controls
    %for injector
    eqs{nComp + 4} = s.rateEq(p,bhp,q,fluid.muL);
    %for producer
    eqs{nComp + 5} = s.ctrlEq(bhp);
    
    
    %%
    % Add source and sink
    for ic = 1 : nComp
       %sink (producer)
       eqs{ic}(wc) = eqs{ic}(wc) - s.compsource(q,sw,F,x,y,L);
    end
    
    % source for total flow equation
    eqs{nComp + 1}(wc) = eqs{nComp + 1}(wc) - q.*(1-sw(wc)).*F(wc);
    
    %sink for total water flow
    eqs{nComp + 2}(wc) = eqs{nComp + 2}(wc) - q.*sw(wc).*epw;
    
end

