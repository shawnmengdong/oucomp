%%
%
% perlOneLineDescription(Returns the Henry coefficients)
% 
%%


function k = getHenryCoef()

   k = [2702.7;  % [He(g)]/[He(l)]  
        2222.22; % [Ne(g)]/[Ne(l)]  
        29.41;   % [CO2(g)]/[CO2(l)]
       ];   

   % Values from wikipedia
   % kHcc = [9.051e-3; 1.101e-2; 0.8317]; % c_aq/c_g (dimensionless).
   % kHpc = = [2702.7; 2222.22; 29.41]; % p/c_aq (in  L*atm/mol).
   
   % switch to SI units
   litre = 1e-3*meter^3;
   k = k*litre*atm;
   
end
