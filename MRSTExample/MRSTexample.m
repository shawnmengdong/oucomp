%% grids and rock
[nx,ny,nz] = deal( 10, 10, 10);
[Lx,Ly,Lz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz ]);
G = computeGeometry(G);
rock = makeRock(G, 30*milli*darcy, 0.3);

%% rock property
cr = 1e-6/barsa;
p_r = 200*barsa;
pv_r = poreVolume(G, rock);
pv = @(p) pv_r .* exp( cr * (p-p_r) );

%% fluid property
mu = 5*centi*poise;
c = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS = 750*kilogram/meter^3;
rho = @(p) rho_r .* exp( c * (p-p_r) );

%% well spec
nperf = 8;
I = repmat(2, [nperf, 1]);
J = (1:nperf).'+1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'producer', 'Dir' , 'x' );

%% initilization
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);

%% discrete operator system
N = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N = N(intInx, :);
n = size(N,1);
C = sparse([(1:n)'; (1:n )'], N,ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x) C*x;
div = @(x) -C'*x;
avg = @(x) 0.5 * (x(N (:,1)) + x(N (:,2)));

hT = computeTrans(G, rock); % Half?transmissibilities
cf = G.cells.faces(:,1);
nf = G.faces.num;
T = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T = T(intInx); % Restricted to interior

gradz = grad(G.cells.centroids(:,3));
v = @(p) -(T/mu).*( grad(p) -g*avg(rho(p)).*gradz );
presEq = @(p,p0,dt) (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0))+ div( avg(rho(p)).*v(p));


%% well model
wc = W(1).cells; % connection grid cells
WI = W(1).WI; % well?indices
dz = W(1).dZ; % depth relative to bottom?hole
p_conn = @(bhp) bhp + g*dz.*rho(bhp); %connection pressures
q_conn = @(p, bhp) WI .* (rho(p(wc)) / mu) .* (p_conn(bhp)-p(wc));
rateEq = @(p, bhp, qS) qS-sum(q_conn(p, bhp))/rhoS;
ctrlEq = @(bhp) bhp-100*barsa;




%% initialize solution and solution spec
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
nc = G.cells.num;
[pIx, bhpIx, qSIx] = deal(1:nc, nc+1, nc+2);

numSteps = 52; % number of time?steps
totTime = 365*day; % total simulation time
dt = totTime / numSteps; % constant time step
tol = 1e-12; % Newton tolerance
maxits = 10; % max number of Newton its

sol = repmat(struct('time', [ ], 'pressure' , [ ], 'bhp', [ ],'qS' , []), [numSteps + 1, 1]);
sol(1) = struct('time', 0, 'pressure' , double(p_ad),'bhp', double(bhp_ad), 'qS', double(qS_ad));

%% time step solving
t = 0; step = 0;
while t < totTime,
    t = t + dt; step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n',step, convertTo(t-dt,day), convertTo(t, day));
    % Newton loop
    resNorm = 1e99;
    p0 = double(p_ad); % Previous step pressure
    nit = 0;
    while (resNorm > tol) && (nit <= maxits)
        %% Newton update
        %calculate residual form of eqs
        eq1 = presEq(p_ad, p0, dt);
        eq1(wc) = eq1(wc)-q_conn(p_ad, bhp_ad);
        eqs = {eq1, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};
        
        
        eq = cat(eqs{:});
        J = eq.jac{1}; % Jacobian
        res = eq.val; % residual(value)
        upd = -(J\res);  %linear solver
        
        % Newton update
        p_ad.val = p_ad.val + upd(1:1000);
        bhp_ad.val = bhp_ad.val + upd(1001);
        qS_ad.val = qS_ad.val + upd(1002);
        
        resNorm = norm(res);
        nit = nit + 1;
        fprintf(' Iteration %3d: Res = %.4e\n', nit, resNorm);
    end
    if nit > maxits
        error('Newton solves did not converge')
    else % store solution
        sol(step+1) = struct('time', t, 'pressure' , double(p_ad),'bhp', double(bhp_ad), 'qS', double(qS_ad));
    end
end

