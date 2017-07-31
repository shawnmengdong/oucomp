%%
%
% perlOneLineDescription(Returns compressibility coefficients)
% 
%%


function omega = omega_l(p, system)

   [nComp, mv, cl, p_ref] = deal(system.nComp, system.fluid.mv, system.fluid.cl, ...
                                 system.fluid.p_ref);
   omega = cell(nComp + 1, 1);
   for ic = 1 : nComp + 1
      omega{ic} = mv*exp(-cl*(p - p_ref));
   end
   
end