%%
% System functions
%  
%  At the beginning of the simulation, we compute the non-dynamic variables such as
%
% * the pore volumes, |pv|,
% * the transmissibilities, |T|, 
% * the divergence and gradient operators, |div| and |grad|,
% * the gravity force contribution, |dz|,
%
%%

function s = sysFunctions(G, rock, W)

   nc = G.cells.num;
   nf = G.faces.num;
   %% 
   % Compute pore volumes
   % 
   s.pv = poreVolume(G, rock);
   s.poro = rock.poro;

   %%
   % Compute the half, and then the full, transmissibilities.
   %
   
   hT = computeTrans(G, rock); %half transmissibilities
   cf = G.cells.faces(:,1); 
   T  = 1 ./ accumarray(cf, 1./hT, [nf, 1]); %full transmissibilities
   s.T = T;

   %% 
   % Set up the discrete divergence operator, |div|. It sums up all signed faces'
   % values in each cell.
   % 
   
   N = double(G.faces.neighbors);
   intInx = all(N ~= 0, 2);
   Nint = N(intInx,:);
   intFace = find(intInx);
   val = reshape(ones(length(intFace),1)*[-1 1],[],1);
    C = sparse([intFace;intFace], Nint(:),val, nf, nc);
    s.grad = @(x) C*x;  %map from cell to face
    s.div = @(x) -C'*x; %map from face to cell
    s.avg = @(x) compute_average(x,G);
    s.faceConcentrations = @(flag,conc_c)faceConcentrations(flag, conc_c, N, intInx,nf, nc);
    
    wc = W(1).cells;
    WI = W(1).WI;
   
    s.q_conn = @(p,bhp,mu) WI./mu.*(bhp-p(wc));  %connection rate
    s.rateEq = @(p,bhp,q,mu) q - sum(s.q_conn(p,bhp,mu));  %rate is the sum of all rate
    s.ctrlEq = @(bhp) bhp -W(1).val;
    
    s.compsource = @(q,sw,F,x,y,L) q.*(1-sw(wc)).*F(wc).*(x(wc).*L(wc)+y(wc).*(1-L(wc))); 
    
end

function conc_f = faceConcentrations(flag, conc_c, N, intInx,nf, nc)
   index        = (1:nf)';
   upCell       = N(:, 2);
   upCell(flag) = N(flag, 1);
   % On the interior cell we use upwind
   M = sparse(index(intInx), upCell(intInx), 1, nf, nc);
   conc_f = M*conc_c;
end

function average= compute_average(x,G)
cellno = gridCellNo(G);  %cell number for each repeated face index
cf = G.cells.faces(:,1); %faces index for each cell, unrepeated
nf_repeat = length(cellno); %number of cell with repeated count
nf = G.faces.num; %number of cell, unrepeated
prop_repeat = zeros(nf_repeat,1);   

for i = 1:nf_repeat
   prop_repeat(i) =  x(cellno(i));
end
accumcount = accumarray(cf,1);%get occurrence of all faces,accumulate the same face
average = accumarray(cf, prop_repeat, [nf, 1])./accumcount; %calculate the average

end