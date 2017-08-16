function G = setupGrid()

%Number of grid cells
nx = 9;
ny = 9;
nz = 4;

%grid dimensions in meter
dx = [0;ones(nx,1)*293.3]*ft;
dy = [0;ones(ny,1)*293.3]*ft;
dz = [0;30;30;50;50]*ft;

%depth at top in meter
z0 = 7315*ft;
x = cumsum(dx);
y = cumsum(dy);
z = cumsum(dz)+z0;

G = tensorGrid(x,y,z);
G = computeGeometry(G);



end