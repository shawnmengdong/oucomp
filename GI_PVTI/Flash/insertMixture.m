function [vapor_y_o,liquid_x_o] = insertMixture(vapor_y,liquid_x,zero_index)

size_o = length(vapor_y)+length(zero_index); %length of output
vapor_y_o =zeros(1,size_o); % initialize output to be all 0
liquid_x_o =zeros(1,size_o); % initialize output to be all 0
index_o = setxor(1:size_o,zero_index); %real index that we kept in checkMixture
vapor_y_o(index_o) = vapor_y; %fill output with data
liquid_x_o(index_o) = liquid_x; %fill output with data

end