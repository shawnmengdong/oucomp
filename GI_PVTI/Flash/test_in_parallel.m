function [success_flag,test_component,z,PT] = test_in_parallel(P_min,P_max,T_min,T_max,options,Comp_Name)
    n = randi(8)+2;
    test_component = randsample(Comp_Name,n);
    z = mynormalize(rand(1,n));
    T = T_min+rand*(T_max-T_min);
    P = P_min+rand*(P_max-P_min);
    [component, ~] = addComponents(test_component);
    thermo = addThermo();
    thermo.EOS = @PREOS;
    mixture = addMixture(component, T, P);
    mixture.mole_fraction = z;
    [success_flag,~,~,~,~,~]=GI_flash(mixture,thermo,options);
    PT = [P,T];

end