function [success_flag,vapor_y,liquid_x,vapor_frac,zc,cubic_time]=vle2fash(mixture,thermo,options)
    eps1 = options.convergence_eps;
    eps2= options.trivial_eps;
    max_itr_RR = options.RRiteration;
    max_itr_outerloop=options.max_outer_loop;
    %Initial estimate for k-values using an empirical equation
    cubic_time = 0;
    ki = wilsonCorrelation(mixture);
    composition = mixture.mole_fraction;
    success_flag =1;
%-------------------------------------------------------------
    [ki_max, max_index] = max(ki);
    % internal step: Correct the maximum ki
    if (ki_max<1)
            ki(max_index) = 1 + ki_max;
    end
    [ki_min, min_index] = min(ki);
    if (ki_min>1)
            ki(min_index) = 0.1;
    end
%--------------------------------------------------------------
    j=0;
    while (1)
        j=j+1;
        % Step 2a: find min and max K-value
        ki_min = min(ki);
        ki_max = max(ki);
        % Step 2b: find min and max vapor fractions
        vap_frac_min = 1/(1-ki_max);
        vap_frac_max = 1/(1-ki_min);
        [vapor_frac, NRflag,~] = RachfordRiceNR(composition, ki, 0.5*(vap_frac_min+vap_frac_max),max_itr_RR);
        % Bi-section
        if ((vapor_frac<vap_frac_min) || (vapor_frac>vap_frac_max) || (NRflag==0))
            f=@(vapor_frac) sum(composition.*(ki-1)./(1+vapor_frac*(ki-1)));
            vapor_frac = bisection(f, vap_frac_min+eps2, vap_frac_max-eps2);
        end  
        
        % Step 4: calculate x and y and normilize
        [liquid_x, vapor_y] = xy_calc(composition, vapor_frac, ki);
        liquid_x = mynormalize(liquid_x);  %!!!!!!!!!!!!!!
        vapor_y = mynormalize(vapor_y); %!!!!!!!!!!!!!!!!
        
        if max(isnan(liquid_x))==1 || max(isnan(vapor_y))==1   %if ki does not converge
            success_flag=3;
            vapor_y = nan;
            liquid_x = nan;
            vapor_frac = nan;
            break;
        end

        [liq_fug, vap_fug,fug_flag,zc] = fugacity(mixture, thermo, liquid_x, vapor_y);
        
        if fug_flag ==0   
            success_flag=4;
            vapor_y = nan;
            liquid_x = nan;
            vapor_frac = nan;
            break;
        end
        
        % Step 6: check the equality of fugacities
        fug_index = (vap_fug~=0);
        error1 = sum(abs(liq_fug(fug_index)./vap_fug(fug_index)-1));
        if (error1<eps1)
            break
        end
        if j>= max_itr_outerloop  %max iteration
            success_flag = 5;
            break;
        end
        
        if (sum(log(ki).^2)<1e-4)  %trival solution
            if (vapor_frac<0)
               vapor_frac=0;
            elseif (vapor_frac>1)
               vapor_frac=1;
            end
            liquid_x = composition;
            vapor_y = composition;
            break
        end

         % Step 7: update K-values
         ki = ki.*liq_fug./vap_fug;
         cubic_time = cubic_time + 2;
    end

end