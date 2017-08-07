function W = setupWells(G,rock)

nperf = 3;
%----producer -----
prod_coord = [2,2];
I_p = repmat(prod_coord(1), [nperf, 1]);
J_p = repmat(prod_coord(2), [nperf, 1]);
K_p = (1:nperf).'+1;
cellInx_p = sub2ind(G.cartDims, I_p, J_p, K_p);
producer_bhp = 4e6;  %pascal
W = addWell([],G,rock,cellInx_p,'type','bhp','val',producer_bhp,'Name','prod1','dir','x');



%--injector -----
% inj_coord = [8,8];
% I_i = repmat(inj_coord(1), [nperf, 1]);
% J_i = repmat(inj_coord(2), [nperf, 1]);
% K_i = (1:nperf).'+1;
% cellInx_p = sub2ind(G.cartDims, I_i, J_i, K_i);
% producer_bhp = 4e6;  %pascal
% W = addWell(W,G,rock,cellInx_p,'type','bhp','val',producer_bhp,'Name','prod1','dir','x');


end