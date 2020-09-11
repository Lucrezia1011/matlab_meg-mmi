% stick all the parameters in an array
% (Note that they are strings)
clear all
close all
clc
%% Compile script
cd ~/matlab/matlab_compiler_test

filename = 'mmi_LTA_powergrid';
% filename = 'mmi_LTA_powergrid_test';

% eval(sprintf('mcc2 -v -m %s.m -R -nojvm singleCompThread ',filename))


%%

% freql=[ 4 8; 8 13; 13 25; 25 40; 40 150]; % do delta separately
freql = [25 40];
latent_vars_name = 'latent_vars.csv';

% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/';
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/';
% gridall = dlmread([data_path,'mni_grid.txt']);
% 
% npoints = '1000';
% 
% nrois =  ceil(nnz(gridall)/str2double(npoints));
% param_list = cell(nrois,1);
% 
% for nn = 1:nrois
%     n = num2str(nn);
%     if size(n,2) == 1
%         n = ['00',n];
%     elseif size(n,2) == 2
%         n = ['0',n];
%     end
%     param_list{nn} = n;
% end


npoints = '266'; % 269 common channels
param_list{1} = '001';

opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
% fit_parameters = X.Properties.VariableNames(3:7);
fit_parameters = X.Properties.VariableNames([4,5,8]);

runcompiled = ['run_',filename,'.sh'];               
compv = 'v96'; % compiler version

cd ~/matlab/matlab_compiler_test
command_list = cell(1,size(freql,1)*length(param_list)*length(fit_parameters));

jj = 0;
% make a command on a new line for each parameter
for ff = 1:size(freql,1)
    freq = sprintf('powersens_%.0f-%.0fHz',freql(ff,1),freql(ff,2));
    meg_data_name = sprintf('%s.txt',freq);

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
                ' test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz '...
                ' && ~/matlab/matlab_compiler_test/%s '...
                ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s;\n'], ... %' mv inds_%s.txt %s%s/lme_%s/;\n'],...
                meg_data_name,data_path,meg_data_name,data_path,latent_vars_name,runcompiled,...
                meg_data_name,latent_vars_name,param_list{ii},npoints,fit_parameter,outpath);
                %param_list{ii},data_path,freq,fit_parameter);
            end
        end
    end
end

command_list(jj+1:end) = [];
command_list = cell2mat(command_list);

% write the commands into a swarm file
file_handle = fopen('mmi_LTA_powersens.swarm','w+');
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
mem = '1';  % gigabytes
threads = '2'; % number of threads
bundles = '1'; % limits number of jobs running at the same time
logfolder = '~/matlab/matlab_compiler_test/swarm_logs';

jobid = evalc(sprintf('!swarm --job-name lmix_powersens --gres lscratch:10 -g %s -t %s -b %s --time 01:00:00 --logdir %s -f mmi_LTA_powersens.swarm --sbatch %s --devel',...
    mem,threads,bundles, logfolder,emailnote))

% try starting swarm from a non interactive session
fprintf('swarm --job-name lmix_powersens --gres lscratch:10 -g %s -t %s -b %s --time 01:00:00 --logdir %s -f mmi_LTA_powersens.swarm --sbatch %s\n',...
    mem,threads,bundles, logfolder,emailnote);

return


%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';

if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat')
end

for ff = 1:5
freq = sprintf('powersens_%.0f-%.0fHz',freql(ff,1),freql(ff,2));

datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/ITI/%s/',freq);
nn =1;

Tp = cell(1,length(fit_parameters));
for ii = 1:length(fit_parameters)
    cd([datapath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
%     Tp{ii} = Xv.tStat;
    Tp{ii} = Xv.pValue;
end

figure(ff); clf; 
set(gcf,'name',freq,'color','w','position',[128  420  1508  402])
for ii = 1:length(fit_parameters)
T = struct;
T.label = sensall;
T.time{1} = 1;
T.sampleinfo = [1 1];


subplot(1,5,ii)
T.trial{1} = Tp{ii}; T.avg = Tp{ii};
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-4 4];
ft_topoplotER(cfg, T)
title(fit_parameters{ii})

end

% saveas(gcf,sprintf('~/matlab/figures/ITI%s.png',freq))
end