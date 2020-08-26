% Made August 4th 2020
% Based on mmi_P300_predict.m
clear all
close all
clc

data_list = [];

meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
Nlist([10,24,26,49,53]) = [];
zz= 0;
for sn = Nlist % all subjects with continuos recordings and latent variables
        
    sdan = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    
    for iiN = 1:3 
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;
            data_list{zz} = data_name;
        end
    end
    
   
end

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% For sensor based analysis
load /data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat
sensall = channels;

%% For grid based analysis
gridres= 5;

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/';
Y = dlmread([data_path,'/meg_trials_evoked_cue.txt']);
Y = reshape(Y, 360,116,size(Y,2));

opts = detectImportOptions([data_path,'latent_vars_evoked_cue.csv']);
Xv = readtable([data_path,'latent_vars_evoked_cue.csv'],opts);

%%
for s = 1:length(data_list)
%% Co-register MRI from fiducial positions

data_name = data_list{s};
sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)
processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
f = data.fsample;


filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


%% Get source localized data on the AAL for this subject
Ym = Y(:,:,Xv.subject == str2double(sub));
Ym = mean(Ym,3);

%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples);

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
if isempty(blockmood_match)
    blockmood_match.sample = [];
    blockmood_match.mood = [];
end
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


bestfit_name = '/data/MBDU/MEG_MMI3/data/behavioral/closed_LTA_coefs.csv';
opts = detectImportOptions(bestfit_name);
bf_pars = readtable(bestfit_name,opts); 
bestfit_sub = bf_pars(bf_pars.Var1 == str2double(sub),:);

g = bestfit_sub.gamma;

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
trials_cue(ttdel) = [];

xi = cue_sample;
mood_cue = Fmood(xi); % interpolated mood timecourse
mood_cue(ttdel) = [];
Scue = repmat(sub,length(mood_cue),1);
gamble_cue = cue_match.choice(cue_match.choice~=0);
gamble_cue(ttdel) = [];

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
sens_cue = squeeze(mean(datas(ind,:,:),2));

outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
save_name = sprintf('%s/cueP300_sens',outpath);


ltvcue = table(Scue,trials_cue',mood_cue',E_cue,Elta_cue,Eltas_cue,gamble_cue',...
    'VariableNames',{'subject','trial','mood','E','E_LTA','E_sum','choice'});
save(save_name,'ltvcue','sens_cue');


%% Co-register MRI

mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

grid =mniLeadfields(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid

%% Beamforming

icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
%     Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*E(1)*eye(size(C)); % 1% max singular value =~ 70*noise


[datay,~]= define_trials(cue_sample, data, tasktime, [-.2,1],0);
ntrialsy = length(datay.trial);
dataym = cell2mat(datay.trial);

L = grid.leadfield(grid.inside);

P = cell(size(L));

for ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    
    % %  G O'Neill method, equivalent to ft
    [v,d] = svd(lf'/Cr*lf);
    d = diag(d);
    jj = 2;
    
    lfo = lf*v(:,jj); % Lead field with selected orientation
    
    w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
    
    ve = w'*dataym;
    ve = reshape(ve,[datay.sampleinfo(1,2),ntrialsy]);  
    
    ver = resample(ve,300,f);
    ver = mean(ver,2);
    c = corr(Ym,ver);
    [~,ind] = max(abs(c));
     
    ve = w'*datavem;
    ve = reshape(ve,[datave.sampleinfo(1,2),ntrials]);  
    
    P{ii} = mean(ve,1)*sign(c(ind));  % P300
    
    if mod(ii,100) == 0
        clc
        fprintf('%s\nSAM running %.1f\n',...
            data_name,ii/length(L)*100)
    end

end

P = cell2mat(P');
Pcue = zeros(size(grid.leadfield,2),size(P,2));
Pcue(grid.inside,:) = P;

%% Save grid data

outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
save_name = sprintf('%s/cueP300_grid',outpath);

ltvcue = table(Scue,trials_cue',mood_cue',E_cue,Elta_cue,Eltas_cue,gamble_cue',...
    'VariableNames',{'subject','trial','mood','E','E_LTA','E_sum','choice'});

save(save_name,'Pcue','ltvcue');
end
return
%% Create combined table for sensor data

ltv_cue = [];

S_cue = [];

r = 0;
for s = 1:length(data_list)
    data_name = data_list{s};
    sub = data_name(5:9);
    outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];

    save_name = sprintf('%s/cueP300_sens',outpath);

    load(save_name)
    r = r+1;
    ltvcue.recording = repmat(r,size(ltvcue,1),1);
    ltv_cue = cat(1,ltv_cue,ltvcue);
 
    S_cue = cat(2,S_cue,sens_cue);   
   
end
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/';

dlmwrite([outpath,'Mcue_P300_30Hzlowpassp.txt'],S_cue)
writetable(ltv_cue,[outpath,'/latent_vars_cuep.csv']);

%% Create combined table for grid data

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/mni_grid.txt');
ltv_cue = [];

P_cue = [];


r = 0;
for s = 1:length(data_list)
    data_name = data_list{s};
    sub = data_name(5:9);
    outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];

    save_name = sprintf('%s/cueP300_grid',outpath);

    load(save_name)
    r = r+1;
    ltvcue.recording = repmat(r,size(ltvcue,1),1);
    ltv_cue = cat(1,ltv_cue,ltvcue);
 
    P_cue = cat(2,P_cue,Pcue);  
   
end
P_cue = P_cue(gridall==1,:);
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/P300/';

dlmwrite([outpath,'BFcue_P300_30Hzlowpassp.txt'],P_cue)
writetable(ltv_cue,[outpath,'/latent_vars_cue.csv']);



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
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
end
sensall = channels;

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/';
freql = {'cue'};
ff = 1;

latent_vars_name = sprintf('latent_vars_%sp.csv',freql{ff});
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
fit_parameters = X.Properties.VariableNames(3:7);

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

T = struct;
T.label = sensall;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position',[176 348 1262 800])

for ii = 1:length(fit_parameters)
    subplot(2,3,ii)
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
saveas(gcf,sprintf('~/matlab/figures/%s.png',freq))

%%
% Plot p-values with multiple comparison correction

figure; set(gcf,'color','w')

for ii = 2:4
    subplot(1,3,ii-1)
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
saveas(gcf,sprintf('~/matlab/figures/%s_pvalue.png',freq))

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