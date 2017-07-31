function [a, b] = simple_mixing_rule(mixture, thermo, ai, bi)
% 

bip = mixture.bip;
mixing_rule_num = thermo.mixingrule;
temperature = mixture.temperature;
x = mixture.mole_fraction;

bipeos = [bip.EOScons]+[bip.EOStdep]*temperature;
N = length(x);
if (mixing_rule_num == 1)  %simple van der Waals mixing
    b = x*bi';
    a=0;
       for i=1:N
         for j=1:N
           a=a+x(i)*x(j)*sqrt(ai(i)*ai(j))*(1-bipeos(i,j));
         end
       end
end

end
    
