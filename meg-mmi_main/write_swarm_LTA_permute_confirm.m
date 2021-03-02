% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_trials_permute_confirm';
% unix(['cp ~/matlab/',filename,'.m ~/matlab/matlab_compiler_test/.'])
% no java virtual machine flag
% eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))

%% Set up parameters
% nrois = 116;

% Hypotheses: 
% Rigth precuneus (68) and paracentral lobule (70) for RPE_LTA and RPE_sum
% Right insula (30) for mood
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

% freq = 'beta';
freq  = 'evoked_outcome';
fit_parameter = 'RPE_LTA'; % RPE_LTA, RPE_sum, mood
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/';
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/';

nroi = 68;

npoints = '360';
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);

fit_parameters = X.Properties.VariableNames(3:end);
% fit_parameters([6,7])= []; % Try for expectation

runcompiled = ['run_',filename,'.sh'];               
compv = 'v98'; % compiler version

out_path = sprintf('%s%s/lme_%s/ROI%d_permute/',data_path,freq,fit_parameter,nroi);
d = dir(out_path);
% 20 permutations per run --> N=100 for 2000 permutations
N = 100 - (length(d)-2);% make a command on a new line for each parameter
% command_list = [];    

randseed = rng('shuffle');

command_list = cell(1,N);
% command_list = cell(1,1*N);


jj = 0;
m = find(strcmp(fit_parameters,fit_parameter));

for nn = 1:N
    % Same random seed for all ROIs
    randseed = rng(randseed.State(2));
    SD = num2str(randseed.Seed);
    jj = jj +1;
        
%      command_list{jj} =  sprintf(['sleep $(echo "${RANDOM}/360" | bc); '...
%         ' export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
%         ' cd /lscratch/$SLURM_JOBID;' ...
%         ' cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/;'...
%         ' test -d /lscratch/$SLURM_JOB_ID/v98 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v98.tar.gz'...
%         ' && ~/matlab/matlab_compiler_test/%s'...
%         ' /lscratch/$SLURM_JOB_ID/v98 %s %s %s %s %s %s;'...
%         ' mv ROI_permute.txt %s%s/lme_%s/ROI%d_permute/%s.txt;\n'],...
%         data_path,meg_data_name,data_path,latent_vars_name,runcompiled,...
%         meg_data_name,latent_vars_name,num2str(nroi),npoints,fit_parameter,SD,...
%         data_path,freq,fit_parameter,nroi,SD);
%     % Add for bundled jobs:  ' rm -rf /lscratch/$SLURM_JOB_ID/*
%     
    command_list{jj} =  sprintf(['sleep $(echo "${RANDOM}/360" | bc); '...
        ' export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
        ' cd /lscratch/$SLURM_JOB_ID; ' ...
        ' if [ -d "%s" ]; then echo "mrc already in scratch"; else cp /usr/local/matlab-compiler/%s.tar.gz /lscratch/$SLURM_JOB_ID/ ' ...
        ' && tar -C /lscratch/$SLURM_JOB_ID -xf %s.tar.gz && rm %s.tar.gz ; fi ; '...
        ' test -f %s || cp %s%s /lscratch/$SLURM_JOB_ID/;'...
        ' test -f %s || cp %s%s /lscratch/$SLURM_JOB_ID/;'...
        ' if [ -f "%s" ]; then echo "y"; else cp ~/matlab/matlab_compiler_test/%s /lscratch/$SLURM_JOB_ID/ '...
        ' && cp ~/matlab/matlab_compiler_test/%s /lscratch/$SLURM_JOB_ID/; fi;',...
        ' /lscratch/$SLURM_JOBID/%s /lscratch/$SLURM_JOB_ID/%s %s %s %s %s %s %s;'...
        ' mv ROI_permute.txt %s/%s.txt;\n'],...
        compv,compv,compv,compv, ...
        meg_data_name,data_path,meg_data_name,...
        latent_vars_name,data_path,latent_vars_name,...
        runcompiled,runcompiled,filename,...
        runcompiled,compv,...
        meg_data_name,latent_vars_name,num2str(nroi),npoints,fit_parameter,SD,...
        out_path,SD);
end



command_list = cell2mat(command_list);

% write the commands into a swarm file
file_handle = fopen(sprintf('mmi_LTA_trials_permute_%s-%d.swarm',fit_parameter,nroi),'w+');
% file_handle = fopen(sprintf('mmi_LTA_trials_%s_missing.swarm',freq),'w+');

fprintf(file_handle,command_list);
fclose(file_handle);


if ~exist(out_path,'dir')
    mkdir(out_path)
end

% Try 1 instance of compiled script
% eval(sprintf('!./%s /usr/local/matlab-compiler/v98 %s %s %s %s %s',runcompiled,meg_data_name,latent_vars_name,'002',npoints,fit_parameter))

%% Run swarm
clc
if ~exist(sprintf('%s/%s',data_path,freq),'dir')
    mkdir(sprintf('%s/%s',data_path,freq))
end
cd ~/matlab/matlab_compiler_test

emailnote = '"--mail-type=FAIL,END"';
% need to include lscratch! see matlab biowulf page
mem = '6';  % gigabytes
threads = '2'; % number of threads
bundles = '5'; %3 % limits number of jobs running at the same time

logfolder = '~/matlab/matlab_compiler_test/swarm_logs';
jobid = evalc(sprintf('!swarm --job-name AAL_permute_%s_%d --gres lscratch:20 -g %s -t %s -b %s --time 4:00:00 --logdir %s -f mmi_LTA_trials_permute_%s-%d.swarm --sbatch %s --devel',...
    fit_parameter,nroi,mem,threads,bundles, logfolder,fit_parameter,nroi,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name AAL_permute_%s_%d --gres lscratch:20 -g %s -t %s -b %s --time 4:00:00 --logdir %s -f mmi_LTA_trials_permute_%s-%d.swarm --sbatch %s\n',...
    fit_parameter,nroi,mem,threads,bundles, logfolder,fit_parameter,nroi,emailnote); % --partition=quick,norm

    