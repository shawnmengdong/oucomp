function s = sysFunctions(system)

s = struct();

G = system.G;
rock = system.rock;
W = system.W;
nc = G.cells.num;
nf = G.faces.num;


%---PoreVolume function--
pv_r = poreVolume(G, rock);
s.pv = @(p) pv_r .* exp(rock.cr * (p - rock.p_r));

%--half,and full transmissibilities---
hT = computeTrans(G, rock); %half transmissibilities
cf = G.cells.faces(:,1); 
T  = 1 ./ accumarray(cf, 1./hT, [nf, 1]); %full transmissibilities
s.T = T;

%--divergence and gradient operator---
%divergence: for a cell,div = face(East)-face(West)+face(North)-face(South)+face(Up)-face(Down)
%gradient: for a face,grad = cell(East)-cell(West) or cell(Up)-cell(Down)or cell(North)-cell(South)
%ignoring boundary flow!
N = double(G.faces.neighbors);
intInx = all(N ~= 0, 2); 
Nint = N(intInx,:);  %interior faces
intFace = find(intInx); %interior faces index
val = reshape(ones(length(intFace),1)*[-1 1],[],1);
C = sparse([intFace;intFace], Nint(:),val, nf, nc);
s.grad = @(x) C*x;  %map from cell to face
s.div = @(x) -C'*x; %map from face to cell
%----------------------------------------------------
s.avg = @(x) compute_average(x,G);
s.faceConcentrations = @(flag,conc_c)faceConcentrations(flag, conc_c, N, intInx,nf, nc);

end


function average= compute_average(x,G)
cellno = gridCellNo(G);  %cell number for each repeated face index, mapping from half face to cell
cf = G.cells.faces(:,1); %faces index for each cell, unrepeated
nf_repeat = length(cellno); %number of faces with repeated count
nf = G.faces.num; %number of cell, unrepeated
prop_repeat = zeros(nf_repeat,1);   

for i = 1:nf_repeat
   prop_repeat(i) =  x(cellno(i));
end
accumcount = accumarray(cf,1);%get occurrence of all faces,accumulate the same face
average = accumarray(cf, prop_repeat, [nf, 1])./accumcount; %calculate the average

end

function conc_f = faceConcentrations(flag, conc_c, N, intInx,nf, nc)
   index        = (1:nf)';
   upCell       = N(:, 1);
   upCell(flag) = N(flag, 2);
   % On the interior cell we use upwind
   M = sparse(index(intInx), upCell(intInx), 1, nf, nc);
   conc_f = M*conc_c;
end