%%
% 
% perlOneLineDescription(Setup the grid from raw data)
%
%%

%%
% load the following data fields:
% * permeability |PermHorix|
% * porosity |Poro|
% * top and bottom depth |Top|, |Bottom|
% * carthesian grid mesh |X|, |Y|

load('RawData');


%%
% plot the top and bottom surface
%

clf
hold on
mesh(X, Y, Bottom, 0*X);
mesh(X, Y, Top, 1+0*X);
hold off


%% 
% Define
%
dims = [99 99 1];
Easting = X*1000;
Northing = Y*1000;
Bottom(Top<Bottom) = Top(Top<Bottom);
Bottom = Bottom';
Top = Top';


physDims=[max(max(Easting))-min(min(Easting)) max(max(Northing))-min(min(Northing)) min(min(Bottom))-max(max(Top))];
g.cartDims = reshape(dims, 1, []);
[X, Y, Z]  = ndgrid(linspace(0, physDims(1), dims(1) + 1), ...
                    linspace(0, physDims(2), dims(2) + 1), ...
                    linspace(0, physDims(3), dims(3) + 1));

Z(:,:,1) = Bottom;
Z(:,:,2) = Top;
actnum = zeros(dims(1:2));
for k=1:2
   actnum = actnum+Z(1:end-1,1:end-1,k)+Z(2:end,1:end-1,k)+Z(1:end-1,2:end,k)+Z(2:end,2:end,k);
end
Z(isnan(Z)) = 500;

lines          = zeros([prod(dims([1, 2]) + 1), 6]);
lines(:,[1,4]) = reshape(X(:,:,[1,end]), [], 2);
lines(:,[2,5]) = reshape(Y(:,:,[1,end]), [], 2);
lines(:,[3,6]) = reshape(Z(:,:,[1,end]), [], 2);

g.COORD = reshape(lines', [], 1);

ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
z   = Z(ind(1), ind(2), ind(3));

g.ZCORN = -z(:); %??
g.ACTNUM =~ isnan(actnum(:));
G = processGRDECL(g, 'RepairZCORN', true);
G = G(1);
G = computeGeometry(G);

rock.perm = interp2(Easting,Northing,PermHoriz*1e-15,G.cells.centroids(:,1),G.cells.centroids(:,2));
rock.perm(rock.perm < 1e-15) = 1e-15;
rock.perm(isnan(rock.perm)) = 1e-15;

rock.poro = interp2(Easting,Northing,Poro,G.cells.centroids(:,1), ...
                  G.cells.centroids(:,2));
rock.poro(rock.poro < 0.2) = 0.2;
rock.poro(isnan(rock.poro)) = 0.2;
