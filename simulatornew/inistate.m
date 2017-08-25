function state0 = inistate(system)
%primary variables: mi (i from 1 to nc), mw , pw
param = system.param;
z0 = param.Datum.z;
p0 = param.Datum.p;
g = param.g;

p0_cow = system.param.Datum.pcow;
p0_cog = system.param.Datum.pcog;
T = system.param.rTemp;
G = system.G;
fluid = system.fluid;
nc = G.cells.num;
z = G.cells.centroids(:,3);

%now assume main pressure is water pressure
p0_w = p0;
p0_o = p0_w+p0_cow;  %po = pw+pcow
p0_g = p0_o+p0_cog;  %pg = po+pcog
comp0 = system.param.mole_fraction;

[~,~,liquid_frac0,molar_density0,rho0,phase_flag0] = flash_calculation(p0_w,T,comp0,system);

delta_z = z-z0;
delta_p = rho0(2)*g*delta_z;
nComp = fluid.Ncomp;

%prelocate our solution:
mi = cell(1,nComp);

pg = p0+delta_p;
po = pg;
pw = pg;

sw = 0.16*ones(nc,1);
sg = 1-sw;
so = zeros(nc,1);
mw = fluid.w.epw*sw;


m_total = sg*molar_density0(2);

for i = 1:nComp
    mi{i} = m_total*comp0(i);
end


%now loop for each cell, find solution
state0.mi = mi;
state0.mw = mw;
state0.pg = pg;
state0.po = po;
state0.pw = pw;
state0.so = so;
state0.sg = sg;
state0.sw = sw;

end