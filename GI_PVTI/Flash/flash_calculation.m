% This function calls flash routine in GI_PVIT

function [liquid_x,vapor_y,liquid_frac,molar_density,rho,phase_flag] = flash_calculation(pressure,Temp,mf,system)
fluid = system.fluid;

%find unique combination of parameters
nComp = fluid.Ncomp;  %number of components
num_cases = size(pressure,1); %number of cases
R = 8.314;
num_phase = 2;
 
%initialize outputs
liquid_x = zeros(num_cases,nComp); %molar fraction in liquid phase
vapor_y = zeros(num_cases,nComp); %molar fraction in vapor phase
molar_density = zeros(num_cases,num_phase); %in mol/m3  
liquid_frac = zeros(num_cases,1); %liquid molar fraction
rho = zeros(num_cases,num_phase);  %in kg/m3
 
 
 
 
 %Assemble total flash matrix
Flash_matrix = zeros(num_cases,nComp+2);   %Flash matrix is [p,t,zi]
Flash_matrix(:,1) = pressure; 
Flash_matrix(:,2) = Temp;
Flash_matrix(:,3:end) = mf;                 

 
%Find the unique configuration
[unique_flash_matrix,ia,ic] = unique(Flash_matrix,'rows');

num_unique = length(ia);  %number of unique configurations
mixture = fluid.mixture;
 
for i = 1:num_unique
    %fluid property calculation
    mixture.pressure = unique_flash_matrix(i,1);
    mixture.temperature = unique_flash_matrix(i,2);
    mixture.mole_fraction = unique_flash_matrix(i,3:end);
    [success_flag,~,y,x,vapor_frac,zc,phase_flag,~,b]=GI_flash(mixture,fluid.thermo,fluid.thermo_opt);
    if success_flag ~=1 
     error('Flash Failed, please debug your flash');
    end

    %assign solution back to outputs
    prev_index = find(ic==i); % original row index in Flash_matrix for the ith row in unique_flash  
    prev_size = length(prev_index);
    liquid_x(prev_index,:) = repmat(x,prev_size,1); 
    vapor_y(prev_index,:) = repmat(y,prev_size,1);
    liquid_frac(prev_index) = 1-vapor_frac; %liquid molar fraction


    mv = zc*R*Temp(i)/pressure(i); %for single configuration molar volume, in m3/mol
    mv(1) = mv(1) - b*[mixture.components.VSE]*x';  %volume shift for liquid

    molar_density(prev_index,:) = repmat(1./mv,prev_size,1);
    rho(prev_index,:) = repmat(1e-3*[mixture.components.MW]*[x',y']./mv,prev_size,1);  %in kg/m3
end

end