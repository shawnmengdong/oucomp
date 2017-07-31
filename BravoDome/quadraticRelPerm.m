%%
%
% perlOneLineDescription(Returns quadratic relative permeability)
% 
%%


function [krL, krG] = quadraticRelPerm(sL)
    krL = sL.*sL;
    krG = (1 - sL).^2;
 end