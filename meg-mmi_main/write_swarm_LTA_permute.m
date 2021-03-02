% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_trials_permute';
% unix(['cp ~/matlab/',filename,'.m ~/matlab/matlab_compiler_test/.'])
% no java virtual machine flag
% eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))

%% Set up parameters
nrois = 116;

% freq = 'beta';
freq  = 'evoked_outcome';
fit_parameter = 'RPE_LTA'; % RPE_LTA, E_sum
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/';
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/';

npoints = '360';
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);

fit_parameters = X.Properties.VariableNames(3:end);
% fit_parameters([6,7])= []; % Try for expectation

runcompiled = ['run_',filename,'.sh'];               
compv = 'v96'; % compiler version

N =10; %1e3; % number of random iterations x 10

% make a command on a new line for each parameter
% command_list = [];    

randseed = rng('shuffle');

command_list = cell(1,N);
% command_list = cell(1,1*N);


jj = 0;
m = find(strcmp(fit_parameters,fit_parameter));

for nn = 1:N
    % Same random seed for all ROIs
    SD = num2str(randseed.Seed);
    randseed = rng(randseed.State(2));
    jj = jj +1;
    
%     command_list{jj} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID; '...
%         '  cd /lscratch/$SLURM_JOBID; if [ -f "%s"] ;  then  echo "data already in lscratch"; ',...
%         '  else cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/; fi ;'...
%         ' test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz'...
%         ' && ~/matlab/matlab_compiler_test/%s'...
%         ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s;'...
%         ' mv ROI_permute2.txt %s%s/lmixmodel_%s/ROI_permute2/%s.txt\n'],...
%         meg_data_name,data_path,meg_data_name,data_path,latent_vars_name,runcompiled,...
%         meg_data_name,latent_vars_name,num2str(nrois),npoints,fit_parameter,SD,...
%         data_path,freq,fit_parameter,SD);
    
     command_list{jj} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
        ' cd /lscratch/$SLURM_JOBID;' ...
        ' cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/;'...
        ' test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz'...
        ' && ~/matlab/matlab_compiler_test/%s'...
        ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s;'...
        ' mv ROI_permute.txt %s%s/lme_%s/ROI_permute/%s.txt;\n'],...
        data_path,meg_data_name,data_path,latent_vars_name,runcompiled,...
        meg_data_name,latent_vars_name,num2str(nrois),npoints,fit_parameter,SD,...
        data_path,freq,fit_parameter,SD);
    % Add for bundled jobs:  ' rm -rf /lscratch/$SLURM_JOB_ID/*
    
end



command_list = cell2mat(command_list);

% write the commands into a swarm file
file_handle = fopen(sprintf('mmi_LTA_trials_permute_%s-%s.swarm',freq,fit_parameter),'w+');
% file_handle = fopen(sprintf('mmi_LTA_trials_%s_missing.swarm',freq),'w+');

fprintf(file_handle,command_list);
fclose(file_handle);

out_path = sprintf('%s%s/lme_%s/ROI_permute/',data_path,freq,fit_parameter);
% out_path = sprintf('%s%s/lmixmodel_%s/ROI_permute2/',data_path,freq,fit_parameter);

if ~exist(out_path,'dir')
    mkdir(out_path)
end

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
mem = '1';  % gigabytes
threads = '1'; % number of threads
bundles = '5'; %3 % limits number of jobs running at the same time

logfolder = '~/matlab/matlab_compiler_test/swarm_logs';
jobid = evalc(sprintf('!swarm --job-name %s_%s --gres lscratch:10 -g %s -t %s -b %s --time 24:00:00 --logdir %s -f mmi_LTA_trials_permute_%s-%s.swarm --sbatch %s --devel',...
    freq,fit_parameter,mem,threads,bundles, logfolder,freq,fit_parameter,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name %s_%s --gres lscratch:10 -g %s -t %s -b %s --time 48:00:00 --logdir %s -f mmi_LTA_trials_permute_%s-%s.swarm --sbatch %s\n',...
    freq,fit_parameter,mem,threads,bundles, logfolder,freq,fit_parameter,emailnote);

    