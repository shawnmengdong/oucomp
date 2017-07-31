% This function calls flash routine in GI_PVIT

function [liquid_x,vapor_y,liquid_frac,molar_density,rho,phase_flag] = flash_calculation(pressure,Temp,z,system)
%find unique combination of parameters
 nComp = system.nComp;  %number of components
 nc = size(pressure,1); %number of cells
 R = 8.314;
 mmC = system.fluid.mmC; %molar mass of components
 num_phase = 2;
 %initialize outputs
 liquid_x = zeros(nc,nComp); %molar fraction in liquid phase
 vapor_y = zeros(nc,nComp); %molar fraction in vapor phase
 molar_density = zeros(nc,num_phase); %in mol/m3  
 liquid_frac = zeros(nc,1); %liquid molar fraction
 rho = zeros(nc,num_phase);  %in kg/m3
 
 
 %Assemble total flash matrix
 Flash_matrix = zeros(nc,nComp+2);   %Flash matrix is [p,t,zi]
 Flash_matrix(:,1) = pressure; 
 Flash_matrix(:,2) = Temp;
 for ic = 1 : nComp
    Flash_matrix(:,ic+2) = z{ic};                 
 end
 
 %Find the unique configuration
 [unique_flash_matrix,ia,ic] = unique(Flash_matrix,'rows');
 
 num_unique = length(ia);  %number of unique configurations
 mixture = system.mixture;
 
 for i = 1:num_unique
     %fluid property calculation
     mixture.pressure = unique_flash_matrix(i,1);
     mixture.temperature = unique_flash_matrix(i,2);
     mixture.mole_fraction = unique_flash_matrix(i,3:end);
     [success_flag,~,y,x,vapor_frac,zc,phase_flag,~]=GI_flash(mixture,system.thermo,system.flashopt);
     if success_flag ~=1 
         error('Flash Failed, please debug your flash');
     end
     
    %assign solution back to outputs
    prev_index = find(ic==i); % original row index in Flash_matrix for the ith row in unique_flash  
    prev_size = length(prev_index);
    liquid_x(prev_index,:) = repmat(x,prev_size,1); 
    vapor_y(prev_index,:) = repmat(y,prev_size,1);
    liquid_frac(prev_index) = 1-vapor_frac; %liquid molar fraction
    molar_rho = pressure(i)/R/Temp(i)./zc; %for single configuration, in mol/m3
    molar_density(prev_index,:) = repmat(molar_rho,prev_size,1);
    rho(prev_index,:) = repmat(mmC*[x',y'].*molar_rho,prev_size,1);  %in kg/m3
 end

end