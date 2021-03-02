% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
clc
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_trials_grid';

% unix(['cp ~/matlab/',filename,'.m ~/matlab/matlab_compiler_test/.'])
eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))

%% Set up parameters
nrois = 10;

param_list = cell(nrois,1);
for nn = 1:nrois
    
    param_list{nn} = num2str(nn);
end

% nrois = 2;
% param_list = cell(nrois,1);
% param_list{1} = '018';
% param_list{2} = '070';

%%

% freq  = 'theta_choice';
% freq = 'tfs_cue' ;
freq = 'evoked_outcome';


latent_vars_name = ['latent_vars_',freq,'.csv'];

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/';

opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);

% X.RPE_abs = abs(X.RPE);
% writetable(X,[data_path,latent_vars_name]);

fit_parameters = X.Properties.VariableNames([3,5,7]);

runcompiled = ['run_',filename,'.sh'];               
compv = 'v96'; % compiler version

cd ~/matlab/matlab_compiler_test
command_list = cell(1,length(param_list)*length(fit_parameters));

jj = 0;
% make a command on a new line for each parameter
for m = 1:length(fit_parameters)

    fit_parameter = fit_parameters{m};
    outpath = [data_path,'/lme_',fit_parameter,'/'];
    
    if ~exist(outpath,'dir')
        mkdir(outpath) 
    end
    
    for ii = 1:length(param_list)
        meg_data_name = ['meg_trials_',freq,'_',param_list{ii},'.txt'];

        jj = jj+1;
        command_list{jj} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
            '  cd /lscratch/$SLURM_JOBID; if [ -f "%s"] ;  then  echo "data already in lscratch"; ',...
            '  else cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/; fi ;'...
            ' test -d /lscratch/$SLURM_JOB_ID/v98 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v98.tar.gz '...
            ' && ~/matlab/matlab_compiler_test/%s '...
            ' /lscratch/$SLURM_JOB_ID/v98 %s %s %s %s %s \n'],...
            meg_data_name,data_path,meg_data_name,data_path,latent_vars_name,runcompiled,...
            meg_data_name,latent_vars_name,param_list{ii},fit_parameter,outpath);
    end
end
command_list(jj+1:end) = [];
command_list = cell2mat(command_list);
% write the commands into a swarm file
file_handle = fopen(sprintf('mmi_LTA_trials_%s_%d.swarm',freq,nrois),'w+');
% file_handle = fopen(sprintf('mmi_LTA_trials_%s_missing.swarm',freq),'w+');

fprintf(file_handle,command_list);
fclose(file_handle);

%     jobid = evalc(sprintf('!swarm --job-name mmiMEGlmix_%s -g 4 --time 01:00:00 -f mmi_LTA_trials_%s.swarm --devel',fit_parameter,fit_parameter))


% Try 1 instance of compiled script
% eval(sprintf('!./%s /usr/local/matlab-compiler/v96 %s %s %s %s %s',runcompiled,meg_data_name,latent_vars_name,'002',npoints,fit_parameter))

%% Run swarm
clc
if ~exist(sprintf('%s/%s',data_path,freq),'dir')
    mkdir(sprintf('%s/%s',data_path,freq))
end
cd ~/matlab/matlab_compiler_test

emailnote = '"--mail-type=FAIL,END"';
% need to include lscratch! see matlab biowulf page
mem = '2';  % gigabytes
threads = '2'; % number of threads
bundles = '1'; % limits number of jobs running at the same time
logfolder = '~/matlab/matlab_compiler_test/swarm_logs';

jobid = evalc(sprintf('!swarm --job-name lmix%s_%d --gres lscratch:10 -g %s -t %s -b %s --time 0:30:00 --logdir %s -f mmi_LTA_trials_%s_%d.swarm --sbatch %s --devel',...
    freq,nrois,mem,threads,bundles, logfolder,freq,nrois,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name lmix%s_%d --gres lscratch:10 -g %s -t %s -b %s --time 1:00:00 --logdir %s -f mmi_LTA_trials_%s_%d.swarm --sbatch %s\n',...
    freq,nrois,mem,threads,bundles, logfolder,freq,nrois,emailnote);
