%%
%
% perlOneLineDescription(Returns quadratic relative permeability)
% 
%%


function [krL, krG,krW] = quadraticRelPerm(sw)
    krL = sw.*sw;
    krW = krL;
    krG = (1 - sw).^2;
 end