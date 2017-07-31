function [vf, RRflag,count] = RachfordRiceNR(composition, ki, vapor_frac_est,max_itr_RR)
% 

    vapor_frac = vapor_frac_est;
    RRflag = 0;
    count = 0;
    while(1)
        count = count + 1;
        f=0;
        dfdv=0;
        N = length(composition);
        for i=1:N
          if (composition(i)~=0)
              f=f+composition(i)*(ki(i)-1)/(1+vapor_frac*(ki(i)-1));
              dfdv=dfdv-composition(i)*(ki(i)-1)^2/(1+vapor_frac*(ki(i)-1))^2;
          end
        end
        dvf = f/dfdv;
        vapor_frac=vapor_frac-dvf;
        if (abs(dvf)<1e-10)
            RRflag = 1;
            break
        end
        if (count>max_itr_RR)
            RRflag = 0;
            break
        end
    end
    vf = vapor_frac;