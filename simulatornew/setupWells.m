function W = setupWells(G,rock)

%Injection well
coord = [1,1];%coordinates for i and j
K_p = [1:2]'; %perforation k coordinates
num_perf = length(K_p);
I_p = repmat(coord(1),[num_perf,1]);
J_p = repmat(coord(2),[num_perf,1]);
Zr = 7330*ft; %reference depth for BHP
cellInx_p = sub2ind(G.cartDims, I_p, J_p, K_p);
W = addWell([],G,rock,cellInx_p,'Type','bhp','val',4000*psia,'Name','I','dir','x','refDepth',Zr);


coord = [7,7];%coordinates for i and j
K_p = [3:4]'; %perforation k coordinates
num_perf = length(K_p);
I_p = repmat(coord(1),[num_perf,1]);
J_p = repmat(coord(2),[num_perf,1]);
Zr = 7400*ft; %reference depth for BHP
cellInx_p = sub2ind(G.cartDims, I_p, J_p, K_p);
W = addWell(W,G,rock,cellInx_p,'Type','bhp','val',500*psia,'Name','P','dir','x','refDepth',Zr);







end