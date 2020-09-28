clear all
close all
clc
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn = 1:16 %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 1.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') ...
                && exist([data_path,name_list(ii).name,'/beamforming/ICA_artifacts.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)   
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
addpath ~/ppyll1/matlab/svdandpca


load /data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat
%%
for s = 1:length(param_list)
%% Co-register MRI from fiducial positions

data_name = param_list{s};
sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

%% Read events

[bv_match,bv] = match_triggers_fc(data_name);

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;


indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;
[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

% xi = gamble_sample;
% mood = Fmood(xi); % interpolated mood timecourse
%   
% trials =  choice_match.bv_index(choice_match.gamble==1)-12;

% Standard model
A = bv.outcomeAmount; % Outcome
A(isnan(A)) = [];
ntrials = length(A);

Es = (bv.winAmount + bv.loseAmount )/2;
Es(isnan(Es)) = [];
RPEs = A - Es;

% LTA model
EltaH = cumsum(A)./(1:ntrials)'; % Expectation, defined by Hanna

RltaH = A - EltaH; % Assume RPE of first trial is 0

g = 0.8;

E_LTA = zeros(ntrials,1);
RPE = zeros(ntrials,1);
for t = 1:ntrials
    E_LTA(t) = sum( g.^(0:(t-1))' .* EltaH(t:-1:1) );
    RPE(t) = sum( g.^(0:(t-1))' .* RltaH(t:-1:1) );
end

% E_LTA = E_LTA(trials);
% RPE = RPE(trials);

% refer to previous trial
E_LTA(2:end+1) = E_LTA;  
E_LTA(1) = (bv.winAmount(13)+bv.loseAmount(13))/2;
E_LTA(end)  =[];
EltaH(2:end+1) = EltaH;
EltaH(1) = (bv.winAmount(13)+bv.loseAmount(13))/2;
EltaH(end) =[];
% E_LTA = E_LTA(trials);
% EltaH = EltaH(trials);
%% Clean data with ICA

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
% cfg.demean = 'yes';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [1 150];
data = ft_preprocessing(cfg);
f = data.fsample;

if exist([processing_folder,'/ICA_artifacts.mat'],'file')
    load([processing_folder,'/ICA_artifacts.mat']);
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data          = ft_rejectcomponent(cfg, comps,data);

%%

filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt

%%
[~,~,ind]= intersect(sensall,data.label);
% include all available trials (even ones when no choice was made)
cue_sample = cue_match.sample(cue_match.choice~=0);
[datave,ttdel]= define_trials(cue_sample, data, tasktime, [.25,.4],0);
ntrials = length(datave.trial);
datavem = cell2mat(datave.trial);

trials =  cue_match.bv_index(cue_match.choice~=0)-12;

% do mood and expectation (all 3 types) influence the P300?
% does not make too much sense to include RPE
E_cue = Es(trials);
E_cue(ttdel) = [];
Elta_cue = EltaH(trials);
Elta_cue(ttdel) = [];
Eltas_cue = E_LTA(trials);
Eltas_cue(ttdel) = [];

trials_cue = trials;

xi = cue_sample;
mood_cue = Fmood(xi); % interpolated mood timecourse
mood_cue(ttdel) = [];
Scue = repmat(sub,length(mood_cue),1);
gamble_cue = cue_match.choice(cue_match.choice~=0);
gamble_cue(ttdel) = [];

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
sens_cue = squeeze(mean(datas(ind,:,:),2));

outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/';
save_name = sprintf('%s%s',outpath,sub);


n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end


ltvcue = table(Scue,trials_cue',mood_cue',E_cue,Elta_cue,Eltas_cue,gamble_cue',...
    'VariableNames',{'subject','trial','mood','E','E_LTA','E_sum','choice'});
save(save_name,'ltvcue','sens_cue');

end

%%
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/';
ltv_cue = [];

S_cue = [];

r = 0;
for s = 1:length(param_list)
    data_name = param_list{s};
    sub = data_name(1:5);
    save_name = sprintf('%s%s',outpath,sub);
    n = str2double(data_name(end-3));
    if ~isnan(n) %check for number at end of filename
        save_name = [save_name,'_',data_name(end-3)];
    else
        save_name = [save_name,'_1'];
    end
    load(save_name)
    r = r+1;
    ltvcue.recording = repmat(r,size(ltvcue,1),1);
    ltv_cue = cat(1,ltv_cue,ltvcue);
 
    S_cue = cat(2,S_cue,sens_cue);   
   
end

dlmwrite([outpath,'Mcue_P300_30Hzlowpassp.txt'],S_cue)
writetable(ltv_cue,[outpath,'/latent_vars_cuep.csv']);

%%
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

cd /data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new
meg = dlmread('meg_trials_evoked_cue.txt');
meg = reshape(meg,[360,116,size(meg,2)]);

opts = detectImportOptions('latent_vars_evoked_cue.csv');
X = readtable('latent_vars_evoked_cue.csv',opts);
time = linspace(-.2,1,360);

figure; 
clf; set(gcf,'color','w')
subplot(221)
plot(time,mean(meg(:,91:end,:),3))
hold on; 
fill([.25 .4 .4 .25],[-1 -1 1 1]*7e-12,[0 0 1],'facealpha',0.1,'edgecolor','none')
ylim([-1 1]*7e-12);  xlim([-.2 1]);
xlabel('time (s)'); ylabel('average evoked response (T)')
title('cerebellum')
subplot(222)
plot(time,mean(meg(:,43:54,:),3))
hold on; 
fill([.25 .4 .4 .25],[-1 -1 1 1]*7e-12,[0 0 1],'facealpha',0.1,'edgecolor','none')
ylim([-1 1]*7e-12); xlim([-.2 1]);
xlabel('time (s)'); ylabel('average evoked response (T)')
title('occipital cortex')
subplot(223)
plot(time,mean(meg(:,59:70,:),3)); 
hold on; 
fill([.25 .4 .4 .25],[-1 -1 1 1]*7e-12,[0 0 1],'facealpha',0.1,'edgecolor','none')
ylim([-1 1]*7e-12); xlim([-.2 1]);
xlabel('time (s)'); ylabel('average evoked response (T)')
title('parietal cortex')
subplot(224)
plot(time,mean(meg(:,79:90,:),3)); 
hold on; 
fill([.25 .4 .4 .25],[-1 -1 1 1]*7e-12,[0 0 1],'facealpha',0.1,'edgecolor','none')
ylim([-1 1]*7e-12);  xlim([-.2 1]);
xlabel('time (s)'); ylabel('average evoked response (T)')
title('temporal cortex')

%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';

if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat')
end

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/';
freql = {'cue'};
ff = 1;

latent_vars_name = sprintf('latent_vars_%sp.csv',freql{ff});
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
fit_parameters = X.Properties.VariableNames(5:6);

freq = sprintf('M%s_P300_30Hzlowpassp',freql{ff});
meg_data_name = sprintf('%s.txt',freq);
meg = dlmread([data_path,meg_data_name]);

outpath = sprintf('%s%s/',data_path,freq);
nn =1;

% including trials where no choice was made lowers significance
Tfit = cell(1,length(fit_parameters));
pfit = cell(1,length(fit_parameters));
for ii = 1:length(fit_parameters)
    cd([outpath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    Tfit{ii} = Xv.tStat;
    pfit{ii} = Xv.pValue;
end
ii=1;
[p,ind] = sort(pfit{ii});
% Significant sensors: 1 cluster on MRT
ss = ind(p'<0.01./(length(sensall):-1:1));

megs = meg(ss,:);
figure; errorbar(mean(megs),std(megs),'.')

T = struct;
T.label = sensall;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position',[176 348 1262 385])

for ii = 1:2
    subplot(1,2,ii)
    T.trial{1} = Tfit{ii}; T.avg = Tfit{ii};
    cfg = [];
    cfg.channel = sensall;
    cfg.layout = 'CTF275_helmet.mat';
    cfg.zlim = [-4 4];
    ft_topoplotER(cfg, T)

    titlename = fit_parameters{ii};
    k = strfind(titlename,'_');
    titlename(k) = ' ';
    title(titlename)
end
%%
% Plot p-values with multiple comparison correction

figure; set(gcf,'color','w')

for ii = 3:4
    subplot(1,2,ii-2)
    p = sort(pfit{ii});
    semilogy(p)
    hold on
    semilogy(0.05./(length(sensall):-1:1))
%     grid on
    xlabel('sensors')
    ylabel('p-value')
    legend('p-values','FDR','location','best')
    N = nnz(p'<0.05./(length(sensall):-1:1));
    
    title(sprintf('%s: p-value of %.0f sensors < 0.05 (FDR)',fit_parameters{ii},N))
end
% saveas(gcf,sprintf('~/matlab/figures/%s_pvalue.png',freq))

megs = megs(1,:);
s = unique(X.subject);
for sn = 1:length(s)
    rr = X.subject == s(sn,:);
    trials = X.trial(rr);
    
    meg_Elta = megs(rr);

    EltaH = cumsum(meg_Elta)./trials'; 
    ntrials = length(trials);
    g = 0.8;
    E_LTA = zeros(ntrials,1);
    for t = 1:ntrials
        E_LTA(t) = sum( g.^(0:(t-1)) .* EltaH(t:-1:1) );
    end

end


%% Plot Linear mixed effects model for grid

mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

freql = {'cue';'choice'}; 
ff = 1;
latent_vars_name = sprintf('latent_vars_%s.csv',freql{ff});
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
fit_parameters = X.Properties.VariableNames(3:7);

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/mni_grid.txt');
T(1:length(fit_parameters)) = {zeros(size(gridall))};
pV(1:length(fit_parameters)) = {zeros(size(gridall))};

param_list = cell(1,15);
for nn = 1:15
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/BF%s_P300_30Hzlowpass/',freql{ff});
    
for ii = 1:length(fit_parameters)
    cd([datapath,'lme_',fit_parameters{ii}])
    X = zeros(nnz(gridall),1);
    pv = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
        pv((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
    end
    T{ii}(gridall==1) = X;
    pV{ii} = pv;
end

%% Plot pValue
freq = sprintf('BF%s_P300_30Hzlowpass',freql{ff});
% meg_data_name = sprintf('%s.txt',freq);
% meg = dlmread([data_path,meg_data_name]);
% 
% M = size(meg,1);
% 
% C = corr(meg');
% lambda = eig(C);
% % Effective number of independent variables
% Meff = 1 + (M-1)*(1 - var(lambda)/M);
% alpha = 1 - (1 - 0.05)^(1/Meff);


figure; set(gcf,'color','w')

for ii = [1,3,4]
    if ii ==1
        subplot(1,3,1)
    else
        subplot(1,3,ii-1)
    end
    p = sort(pV{ii});
    semilogy(p)
    hold on
    semilogy(0.05./(nnz(gridall):-1:1))
%     grid on
    xlabel('voxels')
    ylabel('p-value')
    legend('p-values','FDR','location','best')
    N = nnz(p'<0.05./(nnz(gridall):-1:1));
    
    title(sprintf('%s: p-value of %.0f sensors < 0.05 (FDR)',fit_parameters{ii},N))
end
%% Plot grid
freq = 1;
for ii = [1,3,4]%1:length(fit_parameters)
    sourceant =[];
    
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    sourceant.pow = T{ii};
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'mni';
    
    crang = [-4 4];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
    if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
        cfg.location   = 'max';
    else
        cfg.location   = 'min';
    end
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
    
    ft_sourceplot(cfg, sourceout_Int);
    
    titlename = fit_parameters{ii};
    k = strfind(titlename,'_');
    titlename(k) = ' ';
    
    title(sprintf('%s \npeak t-value %.1f',...
        titlename,max(abs(sourceant.pow(:)))))
    if strcmp(cfg.method,'slice')
        saveas(gcf,sprintf('~/matlab/figures/BF_P300_30Hzlowpass_%s.png',fit_parameters{ii}))
    end
end