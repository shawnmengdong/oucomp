function thermo_model =  addThermo()
% This function creates a thermodynamic model structure for flash

thermo_model = struct();
thermo_model.mixingrule = 1; %van der Waals mixing rule
%thermo_model.activity_model = @NRTL;
thermo_model.EOS = @PREOS;
thermo_model.phase = 1;
thermo_model.fugacity_switch = 1;