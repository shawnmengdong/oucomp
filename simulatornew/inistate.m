function state0 = inistate(system)
%primary variables: mi (i from 1 to nc), mw , pw
z0 = system.param.Datum.z;
p0 = system.param.Datum.p;
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

[~,~,liquid_frac,molar_density,rho,phase_flag] = flash_calculation(p0_w,T,comp0,system);



nComp = fluid.nc;

%prelocate our solution:
mi = cell(1,nComp);
for i = 1:nComp
    mi{i} = zeros(nc,1);
end
mw = zeros(nc,1);
pw = zeros(nc,1);
po = zeros(nc,1);
pg = zeros(nc,1);




%now loop for each cell, find solution
for i = 1:nc
    z1 = z(i);
    





end