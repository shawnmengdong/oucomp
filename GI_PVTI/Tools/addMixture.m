function mixture =  addMixture()
% This function creates a mixture structure
mixture = struct();
component = struct([]);
CNAMES = {'CO2','N2', 'C1', 'C2', 'C3', 'C4-6', 'C7+1', 'C7+2', 'C7+3'};
Tc=5/9*[548.46000 227.16000 343.08000 549.77400 665.64000 806.54054 838.11282 1058.03863 1291.89071]; %R to K
Pc =6894.76*[1071.33111 492.31265    667.78170    708.34238    618.69739 514.92549 410.74956 247.56341 160.41589]; %in psia
Zc =[ .27408 .29115 .28473 .28463 .27748 .27640 .26120 .22706 .20137];
acentric_factor = [ .27408 .29115 .28473 .28463 .27748 .27640 .26120 .22706 .20137];
MW =[44.01000 28.01300 16.04300 30.07000 44.09700 66.86942 107.77943 198.56203 335.19790];
OMEGAA = [.4572355 .4572355 .5340210 .4572355 .4572355 .4572355 .6373344 .6373344 .6373344];
OMEGAB = [.0777961 .0777961 .0777961 .0777961 .0777961 .0777961 .0872878 .0872878 .0872878];
BIP = [0     0      0       0    0    0   0   0  0;
    .0200    0      0       0    0    0   0   0  0;
     .1000  .0360   0       0    0    0   0   0  0;
     .1300  .0500   0       0    0    0   0   0  0;
     .1350  .0800   0       0    0    0   0   0  0;
     .1277  .1002  .092810  0    0    0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0;
     .1000  .1000  .130663  .006 .006 0   0   0  0];
SSHIFT = [0 0 0 0 0 0 0 0 0]; % volume shift parameter
 
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
mixture.components = component;
mixture.bip.EOScons = BIP;   %ESO constant bip
mixture.bip.EOStdep = zeros(N);
end



