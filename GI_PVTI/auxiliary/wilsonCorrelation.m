function K = wilsonCorrelation(mixture)


p = mixture.pressure;
T = mixture.temperature;
critical_pres = [mixture.components.Pc]; % [pa]
critical_temp = [mixture.components.Tc]; %[K]
acentric_fact = [mixture.components.acentric_factor]; %[-]
K = critical_pres/p.*exp(5.37*(1+acentric_fact).*(1-critical_temp/T));
