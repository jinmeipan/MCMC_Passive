%filename_prepare
%
% input:
% -file: the filename for site,month,model-specific choice

function prep_filename(folder0,file)
fid_fn=fopen([folder0,'/temp/FILENAME.txt'],'w');

fprintf(fid_fn,['RunParamsV2_',file,'.txt\n']);
fprintf(fid_fn,['tb_obs_',file,'.txt\n']);
fprintf(fid_fn,['hyperpar_',file,'.txt\n']);
fprintf(fid_fn,['true_theta.txt\n']);
fprintf(fid_fn,['acceptance_',file,'.txt\n']);
fprintf(fid_fn,['post_tb_',file,'.out\n']);
fprintf(fid_fn,['post_theta_',file,'.out']);
fclose(fid_fn);

clear fid_fn

end



