function [pressure,temperature,Comp] = get_flash_input(state,param)

%now assumes the system pressure is gas pressure
pressure = state.pg;

nc = length(pressure);

temperature = param.rTemp*ones(nc,1);

mi = state.mi;

num_comp = length(mi);

Comp = zeros(nc,num_comp);

for i = 1:num_comp
    Comp(:,i) = mi{i};
end

total_m = repmat(sum(Comp,2),1,num_comp);

Comp = Comp./total_m;


end