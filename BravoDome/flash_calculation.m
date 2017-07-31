%%
%
% perlOneLineDescription(Flash calculation)
%
% Given the total concentration of each component |C| and the pressure |p|,
% compute 
% 
% * the concentrations in the gas phase, |cg|, 
% * the concentrations in the liquid phase, |cl|,
% * the liquid saturation, |s|,
% * the total water concentration, |Cw| and 
% * the water concentration in the liquid, |cw|.
%
% In the case we are considering, the flash equation when both phases are present can be
% reduced to a scalar equation, see <perlAddNoteLink() Bravo tutorial note>.
% 
%%

function [cg, cl, cw, s, Cw] = flash_calculation(C, p, system, init_sL)
   
   
   [nComp, k, R, Temp, vp] = deal(system.nComp, system.k, system.R, system.Temp, system.vp);

   tol = 1e-14;
   
   omega_g = R*Temp*ones(nComp, 1);
   min_pressure = sum(bsxfun(@times, omega_g', C), 2) + vp;


   %%
   % There is a minimum pressure, |min_pressure|, under which no liquid phase
   % exists. The cells where this is the case are indexed with |ind_g_unsat|.
   %
   
   ind_g_unsat = p < min_pressure;
   
   %%
   % There is a maximum pressure above which the liquid phase is unsaturated and no gas
   % phase is present. The cells where this is the case are indexed with |ind_l_unsat|.
   % 
   
   ind_l_unsat = sum(bsxfun(@times, C, k') + vp, 2) <= p;

   nc = size(p, 1);
   
   %%
   % If no initial guess for the liquid saturation is given, we start with |sL=0.5|.
   %
   
   if nargin < 4 
      init_sL = 0.5*ones(nc, 1);
   end
   
   s = zeros(nc, 2);
   cg = zeros(nc, nComp);
   cl = zeros(nc, nComp);
   cw = zeros(nc, 2);
   Cw = zeros(nc, 1);

   %%
   % We set the concentrations in the cells where the gas phase is unsaturated. The
   % concentrations in the liquid phase are set to |NaN|.
   %
   
   if any(ind_g_unsat)
      s(ind_g_unsat, 1)  = 1;
      s(ind_g_unsat, 2)  = 0;
      cg(ind_g_unsat, :) = C(ind_g_unsat, :);
      cl(ind_g_unsat, :) = NaN(nnz(ind_g_unsat), 3);
      Cw(ind_g_unsat)    = zeros(nnz(ind_g_unsat), 1);
      cw(ind_g_unsat, :) = zeros(nnz(ind_g_unsat), 2);
   end

   %%
   % We set the concentrations in the cells where the liquid phase is unsaturated. The
   % concentrations in the gas phase are set to |NaN|. We use
   % perlAddLink(computeWaterComp) to compute the water concentrations.
   %
   
   if any(ind_l_unsat)
      s(ind_l_unsat, 1)  = 0;
      s(ind_l_unsat, 2)  = 1;
      cg(ind_l_unsat, :) = NaN(nnz(ind_l_unsat), 3);
      cl(ind_l_unsat, :) = C(ind_l_unsat, :);
      [Cw(ind_l_unsat), cwg, cwl] = computeWaterComp(p(ind_l_unsat), C(ind_l_unsat, :), cl(ind_l_unsat, ...
                                                        :), s(ind_l_unsat, 2), system);
      cw(ind_l_unsat, 1) = cwg;
      cw(ind_l_unsat, 2) = cwl;
   end
   
   %%
   % We compute the liquid saturation by the solving the scalar saturation equation using
   % Newton iterations. Then the functions perlAddLink(computeComposition) and
   % perlAddLink(computeWaterComp) are used to compute the concentrations in the gas
   % and liquid phases.
   
   ind_sat = ~ind_l_unsat & ~ind_g_unsat;
   if any(ind_sat)
      converged = false;
      seq = @(sL)(sat_eq(C(ind_sat, :), p(ind_sat), vp, sL, k, Temp, R, nComp));
      ssat = init_sL(ind_sat);
      while ~converged
         ssat  = initVariablesADI(ssat);
         eqs = seq(ssat);
         dssat = -eqs.val./diag(eqs.jac{1});
         ssat = ssat.val + dssat;
         ssat = min(max(ssat, 0), 1);
         plot(ssat)
         converged = all(abs(eqs.val) <= tol) | all(abs(dssat) < 1e-14);
      end
      s(ind_sat, 1) = 1 - ssat;
      s(ind_sat, 2) = ssat;
      [cg_sat, cl_sat] = computeComposition(C(ind_sat, :), p(ind_sat), ssat, system);
      [Cwsat, cwgsat, cwlsat] = computeWaterComp(p(ind_sat), C(ind_sat, :), cl_sat, ssat, system);
      Cw(ind_sat) = Cwsat;
      cw(ind_sat, 1) = cwgsat;
      cw(ind_sat, 2) = cwlsat;
      cg(ind_sat, :) = cg_sat;
      cl(ind_sat, :) = cl_sat;
   end
end

%%
% This is the scalar equation for the liquid saturation when both phases are present, see
% <perlAddNoteLink() Bravo tutorial note>
%

function y = sat_eq(C, p, vp, s, k, Temp, R, nComp)
   y= -(p/vp - 1);
   for ic = 1:nComp
      y = y + 1/vp*(R*Temp*C(:,ic)./(1 + (R*Temp/k(ic) - 1).*s));
   end
end

