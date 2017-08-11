function  PVTinitialize()



 current_path = fileparts(mfilename('fullpath'));
 addpath([current_path '/Tools']);
 addpath([current_path '/EOS']);
 addpath([current_path '/Flash']);
 addpath([current_path '/auxiliary']);
 addpath([current_path '/MixingRules']);

 disp('Flash Path Completed...');
 
end
