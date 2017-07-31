%%
%
% perlOneLineDescription(Returns the water vapor pressure)
%
%%


function p = vaporPressure(T)

   p = exp(20.386-5132/T)/760; % Formula taken from wikipedia. 1mmHg = 1/760atm
   
   % Switch to SI units
   p = p*atm;
   
end
