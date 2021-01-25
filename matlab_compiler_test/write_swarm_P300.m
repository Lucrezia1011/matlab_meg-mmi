% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
clc
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_powergrid';

% eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))


%%

dimopt = 'M'; %'BF','M'


if strcmp(dimopt,'BF')
    data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/P300/confirm/';
%     gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/P300/mni_grid.txt');
%     gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/mni_grid.txt');
    gridall = dlmread([data_path,'mni_grid.txt']);
    npoints = '1000';
    
    nrois =  ceil(nnz(gridall)/str2double(npoints));
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
    freql = {'cue'};
    
elseif strcmp(dimopt,'M')
    data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/confirm/';
%     freql = {'cue';'choice'};
    freql = {'cue'};
    npoints = '267'; % 266 common channels
    param_list{1} = '001';
end

latent_vars_name = sprintf('latent_vars_%s.csv',freql{1});

opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
% fit_parameters = X.Properties.VariableNames(3:7);
fit_parameters = X.Properties.VariableNames([5,6]);


runcompiled = ['run_',filename,'.sh'];               
compv = 'v98'; % compiler version, changed to 2020a, v98
compiler_path =  ['/data/liuzzil2/matlabCompiler/',filename];

cd ~/matlab/matlab_compiler_test
command_list = cell(1,length(param_list)*length(fit_parameters));

jj = 0;
% make a command on a new line for each parameter
for ff = 1:size(freql,1)
    freq = sprintf('%s%s_P300_30Hzlowpass',dimopt,freql{ff});
%     freq = sprintf('%s%s_P300_30Hzlowpass_planar',dimopt,freql{ff});
    meg_data_name = sprintf('%s.txt',freq);

    latent_vars_name = sprintf('latent_vars_%s.csv',freql{ff});
   
    for m = 1:length(fit_parameters)
        
        fit_parameter = fit_parameters{m};
        outpath = [data_path,freq,'/lme_',fit_parameter,'/'];
        if ~exist(outpath,'dir')
            mkdir(outpath)
        end
        
        for ii = 1:length(param_list)
            
            if ~exist([data_path,freq,'/lme_',fit_parameter,'/inds_',param_list{ii},'.csv'],'file')
            
            jj = jj+1;
            command_list{jj} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
                ' cd /lscratch/$SLURM_JOBID; if [ -f "%s"] ;  then  echo "data already in lscratch"; ',...
                ' else cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/; fi ;'...
                ' test -d /lscratch/$SLURM_JOB_ID/%s || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/%s.tar.gz '...
                ' && ~/matlab/matlab_compiler_test/%s '...
                ' /lscratch/$SLURM_JOB_ID/%s %s %s %s %s %s %s;\n'],... % ' mv inds_%s.txt %s%s/lme_%s/;\n'],...              
                meg_data_name,data_path,meg_data_name,data_path,latent_vars_name,...
                compv,compv,runcompiled,compv,...
                meg_data_name,latent_vars_name,param_list{ii},npoints,fit_parameter,outpath);
                %param_list{ii},data_path,freq,fit_parameter);
%             command_list{jj} =  sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID;'...
%                 ' cd /lscratch/$SLURM_JOBID; cp %s%s /lscratch/$SLURM_JOB_ID/ && cp %s%s /lscratch/$SLURM_JOB_ID/; '...
%                 ' %s/application/%s %s/MATLAB_Runtime/%s'...
%                 ' %s %s %s %s %s %s;\n'],... % ' mv inds_%s.txt %s%s/lme_%s/;\n'],...              
%                 data_path,meg_data_name,data_path,latent_vars_name,...
%                 compiler_path,runcompiled,compiler_path,compv,...
%                 meg_data_name,latent_vars_name,param_list{ii},npoints,fit_parameter,outpath);
%                 %param_list{ii},data_path,freq,fit_parameter);
            
            end
        end
    end
end

command_list(jj+1:end) = [];
command_list = cell2mat(command_list);

% write the commands into a swarm file
file_handle = fopen('mmi_LTA_P300.swarm','w+');
% file_handle = fopen(sprintf('mmi_LTA_trials_%s_missing.swarm',freq),'w+');

fprintf(file_handle,command_list);
fclose(file_handle);

%     jobid = evalc(sprintf('!swarm --job-name mmiMEGlmix_%s -g 4 --time 01:00:00 -f mmi_LTA_trials_%s.swarm --devel',fit_parameter,fit_parameter))


% Try 1 instance of compiled script
% eval(sprintf('!./%s /usr/local/matlab-compiler/v96 %s %s %s %s %s',runcompiled,meg_data_name,latent_vars_name,'002',npoints,fit_parameter))

%% Run swarm
clc

cd ~/matlab/matlab_compiler_test

emailnote = '"--mail-type=FAIL,END"';
% need to include lscratch! see matlab biowulf page
mem = '2';  % gigabytes
threads = '2'; % number of threads
bundles = '1'; % limits number of jobs running at the same time
logfolder = '~/matlab/matlab_compiler_test/swarm_logs';

jobid = evalc(sprintf('!swarm --job-name lmix_%sP300 --gres lscratch:10 -g %s -t %s -b %s --time 01:00:00 --logdir %s -f mmi_LTA_P300.swarm --sbatch %s --devel',...
    dimopt,mem,threads,bundles, logfolder,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name lmix_%sP300 --gres lscratch:10 -g %s -t %s -b %s --time 02:00:00 --logdir %s -f mmi_LTA_P300.swarm --sbatch %s\n',...
    dimopt,mem,threads,bundles, logfolder,emailnote);



