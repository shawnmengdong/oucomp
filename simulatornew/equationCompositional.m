% Residual equations functions
%
% In this function, the residual equations are assembled and evaluated. By using
% _Automatic Differentiation_ (AD), we compute automatically the derivatives of the
% residual. 
%
%
% The mass conservation equations are discretized in space using a Two Point Flux
% Approximation. The transmissibilities are computed before hand when the system is set up
% in the function |setupSystem|. The equations are solved implicitly in time. 
%
function eqs = equationCompositional(state0, state, dt,system)

    s = system.s; % System functions
    fluid = system.fluid;
    T = s.T; % Transmissibilities.
    nComp = fluid.Ncomp; 
    G = system.G;
    %W = system.W;
    %wc = W(1).cells;
    nPhase_max = fluid.nPhase_max; %number of phases except for water
    % 
    % The main variables are converted to AD class objects by the function |initVariablesADI|.  Since they
    % are initialized as main variables, their Jacobian are identity matrices. By operator
    % overloading, any expression which contains these variables will also be an AD object and have a derivative. In this way, the derivatives of
    % the residuals are automatically obtained.
    
    p    = state.pg;
    mi = state.mi;
    mw = state.mw;
    so = state.so;
    sw = state.sw;
    sg = state.sg;
    X = state.X;  %molar fraction in each phase
    phase_frac = state.phase_frac;
    ep = state.ep;
    rho = state.rho;
    %bhp = state.bhp;
    %q = state.q;

    [p,mw,mi{:}] = initVariablesADI(p,mw,mi{:});
    mi0   = state0.mi;
    mw0  = state0.mw;
    
    % Compute the relative permeabilities. 
    kro = fluid.krofn(so);
    krg = fluid.krgfn(sg);
    krw = fluid.krwfn(sw);
    
    % Set up the gravity term. The molar mass are needed.
    g  = 9.8;
    dz = s.grad(G.cells.centroids(:,3));
    rhow   = fluid.w.rhosc; % returns mass density for water.
    epw = fluid.w.epw; %water molar density
    % Compute the mobilities
    mobComp   = {kro.*ep(:,1).*X{:,1}./fluid.muo,krg.*ep(:,2).*X{:,2}./fluid.mug};
    mobW   = krw.*fluid.w.epw./fluid.w.mu;
    
    % Compute the upstream direction of effective potential for each phase.
    dp_phase = cell(1,nPhase_max);
    up_phase = cell(1,nPhase_max);
    
    for ip = 1:nPhase_max
       dp_phase{ip} = s.grad(p) - g*(s.avg(rho(:,ip)).*dz);
       up_phase{ip} = (double(dp_phase{ip})>=0);
    end
    
    dpW = s.grad(p) - g*(rhow*dz);
    upW  = (double(dpW)>=0);
    
    %%
    % Compute the residual of the molar conservation equation for each component 
    for ic = 1 : nComp
        % The function |s.faceConcentrations| computes concentration on the faces given cell concentrations, using upwind directions.
        flux_o = s.faceConcentrations(up_phase{1},mobComp{1}(:,ic)); %liquid flux:mobL*epL*x
        flux_g = s.faceConcentrations(up_phase{2},mobComp{2}(:,ic));
        eqs{ic} = s.pv(p)/dt.*(mi{ic} - mi0{ic}) - s.div(flux_o.*T.*dp_phase{1}+flux_g.*T.*dp_phase{2});
    end
    
    %%
    % Compute the residual of the mass conservation equation for water.
    flux_w = s.faceConcentrations(upW,mobW);
    eqs{nComp + 1} = s.pv(p)/dt.*(mw - mw0) - s.div(flux_w.*T.*dpW);

    %%
    % Compute the residual for saturation balance
    m_total = 0;
    for ic = 1:nComp
        m_total = m_total+mi{ic};
    end
    eqs{nComp + 2} = m_total.*phase_frac(:,1)./ep(:,1)+m_total.*phase_frac(:,2)./ep(:,2)+ mw./fluid.w.epw-1;

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

