function [liquid_x, vapor_y] = xy_calc(composition, vapor_frac, ki)
%calculates the vapor and liquid compositions 
liquid_x = composition./(1+vapor_frac*(ki-1));
vapor_y=ki.*liquid_x;
    