function [stability_flag_l,stability_flag_g,stability_flag] = stabilityTest(mixture, thermo)
% Michelsen stability test; I have used the algorithm described here:
% https://www.e-education.psu.edu/png520/m17_p7.html

trivial_eps = 1e-4;
convergence_eps = 1e-10;
max_itr = 2000;

% extract EOS function
eosf = thermo.EOS;

% switch on the fugacity calculation in thermo structure
thermo.fugacity_switch=1; % switch on

% initialize pseudo second phases
gasmixture = mixture;
liquidmixture = mixture;
gasthermo = thermo;
gasthermo.phase = 2;
liquidthermo = thermo;
liquidthermo.phase = 1;

% extract the total composition and pressure
composition = mixture.mole_fraction;
p = mixture.pressure;

% --------------------------- FIRST TEST ----------------------------------
% calculate the fugacity of the mixture, assuming it is a liquid
[fug_coef,~]=eosf(mixture, liquidthermo);
mixfug = fug_coef.*composition*p; %[Pa]

% Initial estimate for k-values using an empirical equation
ki = wilsonCorrelation(mixture);

% assign large number to error values to begin the loop
conv_error = 1+convergence_eps;
triv_error = 1+trivial_eps;
j = 0;
while (conv_error>convergence_eps) && (triv_error>trivial_eps) && (j<max_itr)
    j = j+1; % loop counter
    % create a vapor-like second phase
    Yi = composition.*ki;
    SV = sum(Yi);

    % normalize the vapor-like mole fractions
    yi = Yi/SV;

    % calculate the fugacity of the vapor-like phase using the thermo structure
    gasmixture.mole_fraction = yi;
    [fug_coef,~]=eosf(gasmixture, gasthermo);
    gasfug = fug_coef.*yi*p; %[Pa]

    % correct K-values
    Ri = mixfug./gasfug/SV;
    ki = ki.*Ri;

    % calculate the convergence and trivial solution error values
    conv_error = sum((Ri-1).^2);
    triv_error = sum(log(ki(ki>0)).^2);
end

% analyze the first test results
if triv_error <= trivial_eps 
    stability_flag_l = 1; % converged to trivial solution
elseif conv_error <= convergence_eps
    stability_flag_l = 2; % converged
else 
    stability_flag_l = 3; % maximum iteration reached
end

% --------------------------- SECOND TEST ---------------------------------
% calculate the fugacity of the mixture, assuming it is a gas
[fug_coef,~]=eosf(mixture, gasthermo);
mixfug = fug_coef.*composition*p; %[Pa]

% Initial estimate for k-values using an empirical equation
ki = wilsonCorrelation(mixture);

% assign large number to error values to begin the loop
conv_error = 1+convergence_eps;
triv_error = 1+trivial_eps;
j = 0;
while (conv_error>convergence_eps) && (triv_error>trivial_eps) && (j<max_itr)
    j = j+1; % loop counter
    % create a liquid-like second phase
    Xi = composition./ki;
    SL = sum(Xi);

    % normalize the liquid-like mole fractions
    xi = Xi/SL;

    % calculate the fugacity of the liquid-like phase using the thermo structure
    liquidmixture.mole_fraction = xi;
    [fug_coef,~]=eosf(liquidmixture, liquidthermo);
    liquidfug = fug_coef.*xi*p; %[Pa]

    % correct K-values
    Ri = liquidfug./mixfug*SL;
    ki = ki.*Ri;

    % calculate the convergence and trivial solution error values
    conv_error = sum((Ri-1).^2);
    triv_error = sum(log(ki(ki>0)).^2);
end

% analyze the second test results
if triv_error <= trivial_eps 
    stability_flag_g = 1; % converged to trivial solution
elseif conv_error <= convergence_eps
    stability_flag_g = 2; % converged
else 
    stability_flag_g = 3; % maximum iteration reached
end


%======================================================================
%interpretation of final results

if stability_flag_g==1 && stability_flag_l ==1
    stability_flag = 1 ; %which is stable
elseif (stability_flag_g==1 && stability_flag_l==3) || (stability_flag_g==3 && stability_flag_l==1) || (stability_flag_g==3 && stability_flag_l==3)
    stability_flag = 2 ; % which is unknown,since it exceeds the maximum iteration
elseif stability_flag_g==3 || stability_flag_l ==3
    switch stability_flag_g
        case 3
            S = SV;
        case 2
            S = SL;
    end
    if S > 1
        stability_flag = 0; % which is unstable
    else
        stability_flag = 2; % which we are not sure because the other phase exceeded max_iter, it may be either > or < = 1 
    end 
else
    if stability_flag_g ==1 && stability_flag_l==2
         S = SV;
    elseif stability_flag_g ==2 && stability_flag_l==1
        S = SL;
    else
        S = max (SL,SV);
    end
    if S > 1
        stability_flag = 0; % which is unstable
    else
        stability_flag = 1; % which is stable
    end 
end

end
