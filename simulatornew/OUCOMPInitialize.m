function OUCOMPInitialize()

 current_path = fileparts(mfilename('fullpath'));
 cd(current_path);
 slash_index = regexp(current_path,'\');
 parrent_path = current_path(1:slash_index(end)-1);
 run([parrent_path,'\GI_PVTI\PVTinitialize.m']);
 run([parrent_path,'\mrst-2017a\startup.m']);
 
end