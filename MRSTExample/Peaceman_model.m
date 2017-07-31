%Peaceman models
[nx,ny,nz] = deal(20,20,5);
G = computeGeometry( cartGrid([nx,ny,nz], [500 500 25]) );
rock = makeRock(G, 100*milli*darcy, .2);
fluid = initSingleFluid('mu',1*centi*poise,'rho',1014*kilogram/meter^3);
hT = computeTrans(G, rock);

%injector, vertical well
W = verticalWell([], G, rock, 1, 1, 1:nz, 'Type', ' rate ' ,'Val' , 3e3/day(), 'Radius', .12*meter, 'name', 'I ' );

W = addWell(W, G, rock, nx : ny : nx*ny, 'Type', 'bhp', 'Val' , 1.0e5, 'Radius', .12*meter, 'Dir' , 'y' , 'name', 'P' );


gravity reset on;
resSol = initState(G, W, 0);
state = incompTPFA(state, G, hT, fluid, 'wells' , W);

cf = accumarray(getCellNoFaces(G),abs(faceFlux2cellFlux(G, state.flux)));