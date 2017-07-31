function phase_flag = phase_Identify(mixture,thermo)


mixing_rule_num = thermo.mixingrule;

critical_pres = [mixture.components.Pc]; %[Pa]
critical_temp = [mixture.components.Tc]; %[K]
acentric_fact = [mixture.components.acentric_factor]; %[-]

p = mixture.pressure; %[Pa]
T = mixture.temperature;
R=8.314;
s1 = 0.623225;  %Huron Vidal
bi = 0.077796*R*critical_temp./critical_pres;
aci=0.457235*(R*critical_temp).^2 ./critical_pres;
mi = 0.37646+(1.54226-0.26992*acentric_fact).*acentric_fact;
Tr = T./critical_temp;
alfai = 1+mi.*(1-sqrt(Tr));   %//alfai=ai^0.5
alfa = alfai.^2.0;           %//alfa=ai
ai = aci .* alfa;
%Q is the parameter for MHV1 and MHV2 mixing rule
%it depends on the EOS
Q = (mixing_rule_num==3)*[-0.53 0]+(mixing_rule_num==4)*[-0.4347 -0.003654];
[a, b] = simple_mixing_rule(mixture, thermo, ai, bi);
A_coef=a*p/(R*T)^2;
B_coef=b*p/(R*T);
poly_coef = [1 -1+B_coef A_coef-B_coef*(2.0+3*B_coef) -B_coef*(A_coef-B_coef*(1+B_coef))];
z_root = roots(poly_coef);
%=======================================================================
if (sum(imag(z_root)~=0)==0)
    z_l = min(z_root);
    z_h = max(z_root);
    critical_fluid_flag = 0;
else
    z_l = z_root(imag(z_root)==0);
    z_h = z_l;
    critical_fluid_flag = 1;
end
del1=1-2^0.5;
del2=1+2^0.5;
difference=(z_h-z_l)*log((z_l-B_coef)/(z_h-B_coef))-A_coef/B_coef/(del2-del1)*log((z_l+del1*B_coef)/(z_l+del2*B_coef)*(z_h+del2*B_coef)/(z_h+del1*B_coef));

if critical_fluid_flag ==0
    if difference >= 0
        phase_flag = 1; %liquid

    elseif difference <0
        phase_flag = 2; %vapor

    end
else
    phase_flag = 3; %critical fluid
end

end
