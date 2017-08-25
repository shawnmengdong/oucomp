function param = setupParam()
param = struct();

%Fluid sample composition %Assumes constant
param.mole_fraction = [.0121 .0194 .65990 .08690 .05910 .09670 .04745 .01515 .00330];


%temperature and pressure at standard condition
param.sc.Temp = (60+459.67)*Rankine;
param.sc.p = 14.7*psia;

%reservoir temperature
param.rTemp = (200+459.67)*Rankine;  %F to K

%datum depth,pressure and capillary pressure
param.Datum.z = 7500*ft;
param.Datum.p = 3550*psia;
param.Datum.pcow = 0*psia;
param.Datum.pcog = 0*psia;

%Seperator Information  3rd stage is stock tank
sep = struct([]);
sep(1).Temp = (80+459.67)*Rankine;
sep(1).p = 815*psia;
sep(2).Temp = (80+459.67)*Rankine;
sep(2).p = 65*psia;
sep(3).Temp = (60+459.67)*Rankine;   
sep(3).p = 14.7*psia;
param.sep = sep;

%Time step
param.dt = 36.5*day;             % Time step
param.total_time = 3650*day;  % Total time
param.g = 9.8;


end