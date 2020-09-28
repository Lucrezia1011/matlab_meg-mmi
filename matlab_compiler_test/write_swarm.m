% stick all the parameters in an array
% (Note that they are strings)

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn =[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = cell(1);
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)   
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

for ii = 1:length(param_list)
    mmi_LTA_aal_prep(param_list{ii})
end

%%
% make a command on a new line for each parameter
command_list = [];
for ii = 1:length(param_list)
    command_list = [command_list ...
        'export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID; '...
        '  test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz '...
        '  && ~/matlab/matlab_compiler_test/run_mmi_LTA_aal_prep.sh '...
        '/lscratch/$SLURM_JOB_ID/v96 '...
        param_list{ii}...
        '\n'];
end


cd ~/matlab/matlab_compiler_test
% write the commands into a swarm file
file_handle = fopen('mmi_LTA_aal_prep.swarm','w+');
% file_handle = fopen('fish.swarm','w+');

fprintf(file_handle,command_list);
fclose(file_handle);