%%
%
% perlOneLineDescription(Returns the gas constant)
% 
%%

function R = getR()
   R = 0.08205736;   % in L*atm*K^(-1)*mol^(-1)

   % switch to SI units:
   litre = 1e-3*meter^3;
   R = R*litre*atm;
   
end
