%%
%
% perlOneLineDescription(Computation of the gas and liquid concentrations)
%
% Computes
%
% * the gas and liquid concentrations |cg| and |cl| 
%
% for given
%
% * total concentrations |C|,
% * pressure |p|,
% * liquid saturation |sL|.
%
% For more details about the equations used here, we refer to <perlAddNoteLink()
% Bravo tutorial note>.
%
%%

function [cg, cl] = computeComposition(C, p, sL, system)
   [nComp, k, R, Temp] = deal(system.nComp, system.k, system.R, system.Temp);
   if iscell(C)
      for ic = 1 : nComp
         cg{ic} = C{ic}./(1 + (R*Temp/k(ic) - 1)*sL);
         cl{ic} = (R*Temp./k(ic)).*cg{ic};
      end
   else
      cg = C./(1 + sL*(R*Temp./k' - 1));
      cl = bsxfun(@times, cg, (R*Temp./k)');
   end      
end

