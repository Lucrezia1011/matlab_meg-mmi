% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_powergrid_permute';
% unix(['cp ~/matlab/',filename,'.m ~/matlab/matlab_compiler_test/.'])
% no java virtual machine flag
% eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))

%% Set up parameters

freql=[25 40];

latent_vars_name = 'latent_vars.csv';

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/confirm/';
gridname = [data_path,'mni_grid.txt'];
% gridname = [data_path,'mni_grid_multiSpheres.txt'];

gridall = dlmread(gridname);


opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
% fit_parameters = X.Properties.VariableNames(4:5); 
fit_parameter= 'E';

runcompiled = ['run_',filename,'.sh'];
compv = 'v98'; % compiler version

cd ~/matlab/matlab_compiler_test

freq = sprintf('powergrid_%.0f-%.0fHz_mu5',freql(1),freql(2));
% freq = sprintf('powergrid_%.0f-%.0fHz_mu5_multiSpheres',freql(1),freql(2));

meg_data_name = sprintf('%s.txt',freq);

out_path = sprintf('%s%s/lme_%s/grid_permute/',data_path,freq,fit_parameter);
d = dir(out_path);
N = 1000 - (length(d)-2);

% make a command on a new line for each parameter
% command_list = [];    

randseed = rng('shuffle');

command_list = cell(1,N);
% command_list = cell(1,1*N);

jj = 0;
% m = find(strcmp(fit_parameters,fit_parameter));

for nn = 1:N
    % Same random seed for all ROIs
    randseed = rng(randseed.State(2));
    SD = num2str(randseed.Seed);
    
    jj = jj +1;
    
%      command_list{jj} =  sprintf(['sleep $(echo "${RANDOM}/360" | bc); '...
%         ' export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
%         ' cd /lscratch/$SLURM_JOBID; cp %s /lscratch/$SLURM_JOB_ID/ && ' ...
%         ' cp %s%s /lscratch/$SLURM_JOB_ID/ && '...
%         ' cp %s%s /lscratch/$SLURM_JOB_ID/;'...
%         ' test -d /lscratch/$SLURM_JOB_ID/%s || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/%s.tar.gz'...
%         ' && ~/matlab/matlab_compiler_test/%s'...
%         ' /lscratch/$SLURM_JOB_ID/%s %s %s %s %s;'...
%         ' mv grid_permute.txt %s%s/lme_%s/grid_permute/%s.txt;\n'],...
%         gridname,data_path,meg_data_name,...
%         data_path,latent_vars_name,...
%         compv,compv,runcompiled,compv,...
%         meg_data_name,latent_vars_name,fit_parameter,SD,...
%         data_path,freq,fit_parameter,SD);
%     % Add for bundled jobs:  ' rm -rf /lscratch/$SLURM_JOB_ID/*
    
    command_list{jj} =  sprintf(['sleep $(echo "${RANDOM}/360" | bc); '...
        ' export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
        ' cd /lscratch/$SLURM_JOB_ID; ' ...
        ' if [ -d "%s" ]; then echo "mrc already in scratch"; else cp /usr/local/matlab-compiler/%s.tar.gz /lscratch/$SLURM_JOB_ID/ ' ...
        ' && tar -C /lscratch/$SLURM_JOB_ID -xf %s.tar.gz && rm %s.tar.gz ; fi ; '...
        ' test -f %s || cp %s%s /lscratch/$SLURM_JOB_ID/;'...
        ' test -f %s || cp %s%s /lscratch/$SLURM_JOB_ID/;'...
        ' if [ -f "%s" ]; then echo "y"; else cp ~/matlab/matlab_compiler_test/%s /lscratch/$SLURM_JOB_ID/ '...
        ' && cp ~/matlab/matlab_compiler_test/%s /lscratch/$SLURM_JOB_ID/; fi;',...
        ' /lscratch/$SLURM_JOBID/%s /lscratch/$SLURM_JOB_ID/%s %s %s %s %s;'...
        ' mv grid_permute.txt %s/%s.txt;\n'],...
        compv,compv,compv,compv, ...
        meg_data_name,data_path,meg_data_name,...
        latent_vars_name,data_path,latent_vars_name,...
        runcompiled,runcompiled,filename,...
        runcompiled,compv,...
        meg_data_name,latent_vars_name,fit_parameter,SD,...
        out_path,SD);
end


command_list = cell2mat(command_list);

% write the commands into a swarm file
file_handle = fopen(sprintf('mmi_LTA_powergrid_permute_%s.swarm',fit_parameter),'w+');
% file_handle = fopen(sprintf('mmi_LTA_trials_%s_missing.swarm',freq),'w+');

fprintf(file_handle,command_list);
fclose(file_handle);


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
mem = '6';  % gigabytes
threads = '2'; % number of threads
bundles = '2'; %3 % limits number of jobs running at the same time

logfolder = '~/matlab/matlab_compiler_test/swarm_logs';
jobid = evalc(sprintf('!swarm --job-name %s_%s --gres lscratch:20 -g %s -t %s -b %s --time 12:00:00 --logdir %s -f mmi_LTA_powergrid_permute_%s.swarm --sbatch %s --devel',...
    freq,fit_parameter,mem,threads,bundles, logfolder,fit_parameter,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name %s_%s --gres lscratch:20 -g %s -t %s -b %s --time 12:00:00 --logdir %s -f mmi_LTA_powergrid_permute_%s.swarm --sbatch %s\n',...
    freq,fit_parameter,mem,threads,bundles, logfolder,fit_parameter,emailnote);

    