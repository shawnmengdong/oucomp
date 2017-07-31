function [success_flag,stability_flag,vapor_y_o,liquid_x_o,vapor_frac,zc,phase_index,cubic_time]=GI_flash(mixture,thermo,options)


%Check for zero mole fraction
[mixture,zero_index] = checkMixture(mixture); %delete 0 components

%Stability test
[~,~,stability_flag] = stabilityTest(mixture,thermo);
success_flag = 1;
if stability_flag == 2
    success_flag = 2;
    vapor_frac=nan;
    vapor_y = nan; %
    liquid_x = nan;
    zc = nan;
    cubic_time=nan;
    phase_index = nan;
elseif stability_flag == 1 % which is single phase
    phase_flag = phase_Identify(mixture);
    if phase_flag ==1 %liquid single phase
        vapor_frac=0;
        thermo.phase = 1;  %liquid phase
        thermo.fugacity_switch = 0; %no need to calculate fugacity for single phase
        [~,~,zl] = thermo.EOS(mixture,thermo);
        zc = [zl,inf];
        vapor_y = zeros(size(mixture.mole_fraction)); %
        liquid_x = mixture.mole_fraction;
        phase_index = [1,0];
        cubic_time=0;
    else  %vapor single phase
        vapor_frac=1;
        thermo.phase = 2;  %vapor phase
        thermo.fugacity_switch = 0; %no need to calculate fugacity for single phase
        [~,~,zv] = thermo.EOS(mixture,thermo);
        zc = [inf,zv];
        vapor_y = mixture.mole_fraction; %
        liquid_x = zeros(size(mixture.mole_fraction));
        phase_index = [0,1];
        cubic_time=0;
    end
else  % unstable
    [success_flag,vapor_y,liquid_x,vapor_frac,zc,cubic_time]=vle2fash(mixture,thermo,options);
    phase_index = [1,1];
end

%Insert 0 component back for output
[vapor_y_o,liquid_x_o] = insertMixture(vapor_y,liquid_x,zero_index);

end