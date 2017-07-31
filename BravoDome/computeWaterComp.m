%%
%
% perlOneLineDescription(Compute the water concentrations)
%
% Computes 
%
% * the total water concentration |Cw|, 
% * the gas water concentration |cwg| and
% * the liquid water concentration |cwl|
%
% for given
%
% * pressure |p|,
% * total concentrations |C|,
% * liquid concentrations |cl|,
% * liquid saturation |sL|.
%
% For more details about the equations used here, we refer to <perlAddNoteLink()
% Bravo tutorial note>.
%
%%

function [Cw, cwg, cwl] = computeWaterComp(p, C, cl, sL, system)
   [nComp, vp, R, Temp] = deal(system.nComp, system.vp, system.R, system.Temp);
   cwg = vp/(R*Temp);
   omega_liq = omega_l(p, system);
   if iscell(C)
      cwl = omega_liq{1}.*cl{1};
      for ic = 2 : nComp
         cwl = cwl + omega_liq{ic}.*cl{ic};
      end
      cwl = 1./omega_liq{nComp + 1}.*(1-cwl);
      Cw = cwg.*(1 - sL) + cwl.*sL;
   else
      omega_liq = cell2mat(omega_liq);
      cwl = sum(bsxfun(@times, omega_liq(1:nComp)', cl), 2);
      cwl = 1./omega_liq(nComp + 1).*(1-cwl);
      Cw = cwg.*(1 - sL) + cwl.*sL;
   end
end
