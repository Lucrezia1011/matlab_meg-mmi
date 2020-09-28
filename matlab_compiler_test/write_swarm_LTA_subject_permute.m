% stick all the parameters in an array
% (Note that they are strings)

%% Compile script
filename = 'mmi_LTA_subjects_permute';
% unix(['cp ~/matlab/',filename,'.m ~/matlab/matlab_compiler_test/.'])
eval(['mcc2 -m -R singleCompThread ',filename])

%% Set up parameters
nrois = 116;
% nrois = 269;
param_list = cell(nrois,1);
for nn = 1:nrois
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

%%
stype = 'sensors';
% stype = 'aal_mu5max_Z';

freq  = ['evoked_',stype];
npoints = '1200';

        

latent_vars_name = 'latent_vars.csv';

% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/';
data_path = '/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/';
load([data_path,'mmiRN_aal_mu5max_Z.mat'],'Nn','Np','Npn');

opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);

fit_parameters = X.Properties.VariableNames(2:(end-3));

runcompiled = ['run_',filename,'.sh'];               
compv = 'v96'; % compiler version

cd ~/matlab/matlab_compiler_test

N = 1e3; % number of random iterations x 10


meg_data_names = {['positive_',stype,'.txt'];['negative_',stype,'.txt'];...
    ['posneg_',stype,'.txt']};


Nnames = {'Npos';'Nneg';'Nposneg'};
outname = {[freq,'_pos'];[freq,'_neg'];[freq,'_posneg']};
for d = 1:length(meg_data_names)
    
    n = 0;
 
    if ~exist([data_path,meg_data_names{d}],'file')
    % reshape data with correct dimensions
        meg = dlmread([data_path, meg_data_name0]);
        meg = reshape(meg,[nrois,str2double(npoints),size(meg,2)]);
        meg = permute(meg,[2,1,3]);
        meg = reshape(meg,[nrois*str2double(npoints),size(meg,3)]);
        dlmwrite([data_path,meg_data_name],meg)
    end
    
    
    if ~exist(sprintf('%s/%s',data_path,outname{d}),'dir')
        mkdir(sprintf('%s/%s',data_path,outname{d}))
    end
    
    % make a command on a new line for each parameter
    command_list =  cell(1,length(param_list)*N);

    
    for m = 1:length(fit_parameters)
        fit_parameter = fit_parameters{m};
        for ii = 1:length(param_list)
%             fileroi = sprintf('%s/%s/glmmodel_%s/ROI_%s.csv',data_path,outname{d},fit_parameter,param_list{ii});           
%             if ~exist(fileroi,'file')
                n = n +1;
                command_list{n} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID; '...
                    '  test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz '...
                    '  && ~/matlab/matlab_compiler_test/%s '...
                    ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s %s %s\n'],runcompiled,...
                    meg_data_names{d},latent_vars_name,param_list{ii},npoints,Nnames{d},fit_parameter,outname{d},data_path);
                
%             end            
        end    
    end
    
    command_list = cell2mat(command_list);

    file_handle = fopen(sprintf('mmi_LTA_subjects_%s_permute.swarm',outname{d}),'w+');

    fprintf(file_handle,command_list);
    fclose(file_handle);

    
end
% write the commands into a swarm file

%     jobid = evalc(sprintf('!swarm --job-name mmiMEGlmix_%s -g 4 --time 01:00:00 -f mmi_LTA_trials_%s.swarm --devel',fit_parameter,fit_parameter))


% Try 1 instance of compiled script
% eval(sprintf('!./%s /usr/local/matlab-compiler/v96 %s %s %s %s %s',runcompiled,meg_data_name,latent_vars_name,'002',npoints,fit_parameter))

%% Run swarm

cd ~/matlab/matlab_compiler_test

emailnote = '"--mail-type=FAIL,END"';
% need to include lscratch! see matlab biowulf page
mem = '1';  % gigabytes
threads = '2'; % number of threads
bundles = '50'; % limits number of jobs running at the same time
logfolder = '~/matlab/matlab_compiler_test/swarm_logs';
jobid = evalc(sprintf('!swarm --job-name glm%s --gres lscratch:10 -g %s -t %s -b %s --time 00:30:00 --logdir %s -f mmi_LTA_subjects_%s.swarm --sbatch %s',...
    freq,mem,threads,bundles,logfolder,freq,emailnote));


    