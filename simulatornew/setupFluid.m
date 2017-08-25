function fluid = setupFluid()

fluid = struct();
component = struct([]);

%Component names
CNAMES = {'CO2','N2', 'C1', 'C2', 'C3', 'C4-6', 'C7+1', 'C7+2', 'C7+3'};

%Critical temperatures for each component in K
Tc=[548.46000 227.16000 343.08000 549.77400 665.64000 806.54054 838.11282 1058.03863 1291.89071]*Rankine; %R to K

%Critical pressure for each component in Pa
Pc =[1071.33111 492.31265    667.78170    708.34238    618.69739 514.92549 410.74956 247.56341 160.41589]*psia;

%Critical Compressibility factor
Zc =[ .27408 .29115 .28473 .28463 .27748 .27640 .26120 .22706 .20137];

%Accentric Factor
acentric_factor = [ .27408 .29115 .28473 .28463 .27748 .27640 .26120 .22706 .20137];

%Molecular Weight
MW =[44.01000 28.01300 16.04300 30.07000 44.09700 66.86942 107.77943 198.56203 335.19790];

%The coefficient for EOS calculation
OMEGAA = [.4572355 .4572355 .5340210 .4572355 .4572355 .4572355 .6373344 .6373344 .6373344];
OMEGAB = [.0777961 .0777961 .0777961 .0777961 .0777961 .0777961 .0872878 .0872878 .0872878];

%Binary Interaction Parameter
BIP = [0     0      0       0    0    0   0   0  0;
    .0200    0      0       0    0    0   0   0  0;
     .1000  .0360   0       0    0    0   0   0  0;
     .1300  .0500   0       0    0    0   0   0  0;
     .1350  .0800   0       0    0    0   0   0  0;
     .1277  .1002  .092810  0    0    0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0];
 
%Volume Shift Parameter (Se) 
SSHIFT = [0 0 0 0 0 0 0 0 0];

N = length(CNAMES);
for n = 1:N
    component(n).name = CNAMES{n};
    component(n).MW = MW(n);
    component(n).Tc = Tc(n);
    component(n).Pc = Pc(n);
    component(n).Zc = Zc(n);
    component(n).acentric_factor = acentric_factor(n);
    component(n).OMEGAA = OMEGAA(n);
    component(n).OMEGAB = OMEGAB(n);
    component(n).VSE = SSHIFT(n);  %dimensionless shift parameter SE
end

%Water Saturation Functions Table,  Sw ---Krw---Pcow=Po-Pw (Psia)
SWFN =[ 0.16  0      50;
        0.18  0      41;
        0.20  0.002  32;
        0.24  0.010  21;
        0.28  0.020  15.5;
        0.32  0.033  12.0;
        0.36  0.049  9.2;
        0.40  0.066  7.0;
        0.44  0.090  5.3;
        0.48  0.119  4.2;
        0.52  0.150  3.4;
        0.56  0.186  2.7;
        0.60  0.227  2.1;
        0.64  0.277  1.7;
        0.68  0.330  1.3;
        0.72  0.390  1.0;
        0.76  0.462  0.7;
        0.8   0.540  0.5;
        0.84  0.620  0.4;
        0.88  0.710  0.3;
        0.92  0.800  0.2;
        0.96  0.900  0.1;
        1.00  1.000  0.0];
SWFN(:,3) = SWFN(:,3)*psia;  %Convert psia to pa

%Gas Saturation Functions Table,  Sg ---Krg---Pcog=Pg-Po (Psia)
SGFN = [0.00  0.000  0.0;
        0.04  0.005  0.1;
        0.08  0.013  0.2;
        0.12  0.026  0.3;
        0.16  0.040  0.4;
        0.20  0.058  0.5;
        0.24  0.078  0.6;
        0.28  0.100  0.7;
        0.32  0.126  0.8;
        0.36  0.156  0.9;
        0.40  0.187  1.0;
        0.44  0.222  1.1;
        0.48  0.260  1.2;
        0.52  0.300  1.3;
        0.56  0.349  1.4;
        0.60  0.400  1.5;
        0.64  0.450  1.6;
        0.68  0.505  1.7;
        0.72  0.562  1.8;
        0.76  0.620  1.9;
        0.80  0.680  2.0;
        0.84  0.740  2.1];
SGFN(:,3) = SGFN(:,3)*psia;  %Convert psia to pa

%Oil Saturation Functions Table,  So ---Kro(only oil and water)---Kro(gas
%oil and water)
SOF3 = [0.00  0.000  0.000;
        0.04  0.000  0.000;
        0.08  0.000  0.000;
        0.12  0.000  0.000;
        0.16  0.000  0.000;
        0.20  0.000  0.000;
        0.24  0.000  0.000;
        0.28  0.005  0.005;
        0.32  0.012  0.012;
        0.36  0.024  0.024;
        0.40  0.040  0.040;
        0.44  0.060  0.060;
        0.48  0.082  0.082;
        0.52  0.112  0.112;
        0.56  0.150  0.150;
        0.60  0.196  0.196;
        0.68  0.315  0.315;
        0.72  0.400  0.400;
        0.76  0.513  0.513;
        0.80  0.650  0.650;
        0.84  0.800  0.800];
    
    
%function from table
fluid.swfn = @(pcow) interp1(SWFN(:,3),SWFN(:,1),pcow,'pchip',0.16);
fluid.krwfn = @(sw) interp1(SWFN(:,1),SWFN(:,2),sw);
fluid.krgfn = @(sg) interp1(SGFN(:,1),SGFN(:,2),sg,'pchip',0.74);
fluid.krofn = @(so) interp1(SOF3(:,1),SOF3(:,2),so,'pchip',0.8);


%Water property
fluid.w.pr = 3550*psia;  %reference pressure
fluid.w.Bw = 1;  %formation volume factor at reference pressure
fluid.w.cw = 0.000003/psia; %water compressibility
fluid.w.mu = 0.31*centi*poise; %water viscosity at reference pressure
fluid.w.rhosc = 63*pound/(ft^3); %water density at standard condition
fluid.w.MW = 18/1000;  %molecular weight in kg/mol
fluid.w.epw = fluid.w.rhosc/fluid.w.MW;



%Assemble fluid
fluid.mixture.components = component;
fluid.mixture.bip.EOScons = BIP;   %EOS constant bip
fluid.mixture.bip.EOStdep = zeros(N);
fluid.Ncomp = N;

%Thermo and options
fluid.thermo = addThermo();
fluid.thermo.EOS = @PREOS;
fluid.thermo_opt.convergence_eps = 1e-12;   %convergence tolerance for fugacity
fluid.thermo_opt.trivial_eps = 1e-3;     %trivial shift for bisection algorithm
fluid.thermo_opt.RRiteration = 200;   %maximum number of Rachford Rice iteration using Newton's method
fluid.thermo_opt.max_outer_loop = 1000;   %max number of fugacity updates

%maximum number of phases (except water)
fluid.nPhase_max=2;
fluid.muo = 1e-3;  % Liquid
fluid.mug = 1e-5;  % Gas

end