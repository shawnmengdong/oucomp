%%
% 
% perlOneLineDescription(Setup the system)
%  
%  At the beginning of the simulation, we compute the non-dynamic variables such as
%
% * the pore volumes, |pv|,
% * the transmissibilities, |T|, 
% * the discrete differential operators, |div| and |grad|,
% * the gravity force contribution, |dz|,
%
%%

function s = setupSystem(G, rock, bc, param)

   cf = G.cells.faces(:,1);
   nf = G.faces.num;
   nc = G.cells.num;
   
   %% 
   % Compute pore volumes
   % 
   
   s.pv = poreVolume(G, rock);
   s.poro = rock.poro;

   %%
   % Compute the half, and then the full, transmissibilities.
   %
   
   T = computeTrans(G, rock);
   cf = G.cells.faces(:,1);
   nf = G.faces.num;
   T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
   s.T = T;

   %% 
   % Set up the discrete divergence operator, |div|. It sums up all signed faces'
   % values in each cell.
   % 
   
   N = double(G.faces.neighbors);
   index = 1:nf';
   faces1 = N(:, 1) ~= 0;
   faces2 = N(:, 2) ~= 0;
   C1  = sparse(index(faces1), N(faces1, 1), ones(nnz(faces1),1), nf, nc);
   C2  = sparse(index(faces2), N(faces2, 2), ones(nnz(faces2),1), nf, nc);
   C = C1 - C2;
   s.div  = @(x)C'*x;

   %%
   % Set up the discrete gradient operator, |grad|. We compute the differences of cell values
   % across each face. It is a linear mapping from cells' to faces' values.
   
   index = 1:nf;
   interior = prod(N, 2)~=0;

   C1_interior = sparse(index(interior), N(interior, 1), ones(nnz(interior), 1), nf, nc);
   C2_interior = sparse(index(interior), N(interior, 2), ones(nnz(interior), 1), nf, nc);

   
   %%
   % Compute the boundary contribution to the gradient operator. They corresponds to the
   % external faces where Dirichlet conditions are given. We are careful to use the
   % correct signs.
   % 

   is_dirichlet_faces1 = N(bc.dirichlet.faces, 1) ~= 0;
   is_dirichlet_faces2 = ~is_dirichlet_faces1;

   dirichlet_faces1 = bc.dirichlet.faces(is_dirichlet_faces1);
   dirichlet_faces2 = bc.dirichlet.faces(is_dirichlet_faces2);

   C1_exterior = sparse(index(dirichlet_faces1), ...
                        N(dirichlet_faces1, 1), ...
                        ones(numel(dirichlet_faces1), 1), nf, nc);
   C2_exterior = sparse(index(dirichlet_faces2), ...
                        N(dirichlet_faces2, 2), ...
                        ones(numel(dirichlet_faces2), 1), nf, nc);

   %%
   % The gradient operator is the sum of internal and boundary contributions.
   %
   
   C = C1_interior + C1_exterior - (C2_interior + C2_exterior);

   pressure_bc = sparse(nf, 1);
   pressure_bc(dirichlet_faces1) = - bc.dirichlet.pressure(is_dirichlet_faces1);
   pressure_bc(dirichlet_faces2) = + bc.dirichlet.pressure(is_dirichlet_faces2);

   s.p_grad = @(p)(C*p + pressure_bc);  

   s.grad = @(val, bc_val)(grad(val, bc_val, nf, C, ...
                                is_dirichlet_faces1, dirichlet_faces1, ... 
                                is_dirichlet_faces2, dirichlet_faces2));

   %%
   % Set up the gravity term.
   %
   
   z = G.cells.centroids(:, 3);
   fz = G.faces.centroids(:, 3);
   s.dz = s.grad(z, fz);

   %% 
   % Define function to compute concentrations on faces using cells' values and upwind
   % directions. Here, the _nature_ of the boundary faces (Dirichlet or note) are not
   % allowed to change, although the _values_ on these faces may change.
   %
   
   s.faceConcentrations = @(flag, conc_c, bc_conc) ...
       faceConcentrations(flag, conc_c, bc_conc, N, interior, dirichlet_faces2, ...
                          dirichlet_faces1,  bc, nf, nc);
   
   s.N = N;
   s.G = G;

end

function dval = grad(val, bc_val, nf, C, ...
                     is_dirichlet_faces1, dirichlet_faces1, ...
                     is_dirichlet_faces2, dirichlet_faces2)

   signed_bc_val = sparse(nf, 1);
   signed_bc_val(dirichlet_faces1) = - bc_val(is_dirichlet_faces1);
   signed_bc_val(dirichlet_faces2) = + bc_val(is_dirichlet_faces2);
   dval = C*val + signed_bc_val;

end

function conc_f = faceConcentrations(flag, conc_c, bc_conc, N, interior, ...
                                                 dirichlet_faces2, dirichlet_faces1, ...
                                                 bc, nf, nc)
   index        = (1:nf)';
   upCell       = N(:, 2);
   upCell(flag) = N(flag, 1);

   % On the interior cell we use upwind
   Mint = sparse(index(interior), upCell(interior), 1, nf, nc);

   logical_dirichlet_faces1 = zeros(nf, 1);
   logical_dirichlet_faces1(dirichlet_faces1) = 1;
   logical_dirichlet_faces1 = logical(logical_dirichlet_faces1);
   logical_dirichlet_faces2 = zeros(nf, 1);
   logical_dirichlet_faces2(dirichlet_faces2) = 1;
   logical_dirichlet_faces2 = logical(logical_dirichlet_faces2);

   external_faces1 = N(:,2)==0;
   external_faces2 = N(:,1)==0;


   % On the exterior faces where no Dirichlet conditions are given we take the value given in
   % the interior cell.
   Mext1 = sparse(index(external_faces1 & ~logical_dirichlet_faces1), ...
                  N(external_faces1 & ~ logical_dirichlet_faces1, 1), 1, nf, nc);
   Mext2 = sparse(index(external_faces2 & ~logical_dirichlet_faces2), ...
                  N(external_faces2 & ~ logical_dirichlet_faces2, 2), 1, nf, nc);

   % On the Dirichlet boundary cells we use upwind, taking the values from boundary
   % conditions when needed We assume flag is logical

   assert(islogical(flag), 'Upstream indices must be given as logical');
   Mdir1 = sparse(index(flag & logical_dirichlet_faces1), ...
                  N((flag & logical_dirichlet_faces1), 1), 1, nf, nc);
   Mdir2 = sparse(index(~flag & logical_dirichlet_faces2), ...
                  N((~flag & logical_dirichlet_faces2), 2), 1, nf, nc);

   M = Mint + Mext1 + Mext2 + Mdir1 + Mdir2;

   % Values of saturation from Dirichlet boundary conditions.
   dconc_all = sparse(nf, 1);
   dconc_all(bc.dirichlet.faces) = bc_conc; 
   dconc = sparse(nf, 1);
   dconc(flag & logical_dirichlet_faces2)  = dconc_all( flag & logical_dirichlet_faces2);
   dconc(~flag & logical_dirichlet_faces1) = dconc_all(~flag & logical_dirichlet_faces1);
   
   conc_f = M*conc_c + dconc;
   
end