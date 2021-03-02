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

%%

roiopt = 'grid';

gridres= 5;

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/';
Y = dlmread([data_path,'/meg_trials_evoked_cue.txt']);
Y = reshape(Y, 360,116,size(Y,2));


opts = detectImportOptions([data_path,'latent_vars_evoked_cue.csv']);
Xv = readtable([data_path,'latent_vars_evoked_cue.csv'],opts);


%%
load /data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat

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

Ym = Y(:,:,Xv.subject == str2double(sub));
Ym = mean(Ym,3);

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

%% Co-register MRI

mri_name = [sub,'_anat+orig.BRIK'];

if ~exist(mri_name,'file')
    unix(['gunzip ',mri_name])
end

mri = ft_read_mri(mri_name,'dataformat','afni_brik');

tagset_shape = mri.hdr.TAGSET_NUM;
tagset_coord = mri.hdr.TAGSET_FLOATS;
tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa

tagset_p = zeros(1,3);  % Ideal orientation {RL; PA; IS}
for ii =1:3
    if strcmp(mri.hdr.Orientation(ii,:),'AP') || strcmp(mri.hdr.Orientation(ii,:),'PA')
        tagset_p(ii) = 2;
    elseif strcmp(mri.hdr.Orientation(ii,:),'LR') || strcmp(mri.hdr.Orientation(ii,:),'RL')
        tagset_p(ii) = 1;
    elseif strcmp(mri.hdr.Orientation(ii,:),'SI') || strcmp(mri.hdr.Orientation(ii,:),'IS')
        tagset_p(ii) = 3;
    end
end

m = [   -1  0   0   mri.dim(1)
    0   -1  0   mri.dim(2)
    0   0   1   1
    0   0   0   1] ;


tagset_coord = tagset_coord(tagset_p,:)'; % fiducials have shuffled coordinates

mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin

mri.transform = mri.transform/m;
fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mri = ft_volumerealign(cfg,mri);

if ~exist([sub,'_coreg.nii'],'file')
    writebrik([sub,'_coreg'],mri);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Segment MRI
if ~exist([processing_folder,'/headmodel.mat'],'file')
    cfg = [];
    cfg.output  = 'brain';
    segmentmri = ft_volumesegment(cfg,mri);
    
    % Head model
    
    cfg = [];
    cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, segmentmri);
    
    save([processing_folder,'/headmodel.mat'],'vol')
else
    load([processing_folder,'/headmodel.mat']);
end
sens = ft_read_sens(data_name,'senstype','meg');

%% AAL atlas
%     gridres = 5; % resolution of beamformer grid in mm

% Load fieldtrip 10mm MNI grid
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
template_grid = sourcemodel;
atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,sourcemodel.unit);

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

clear sourcemodel

%% Sourcemodel warp MNI grid

% sourcemodel based on 5mm grid MNI brain
cfg = [];
cfg.mri = mri;
cfg.warpmni = 'yes';
cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
cfg.unit      = 'm';
cfg.nonlinear = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);
locs = sourcemodel.pos;
if  ~strcmp(roiopt,'grid')
    % Find location of AAL ROIs
    R = length(sourcemodelAAL.tissuelabel);
    locs = zeros(R,3);
    locsAAL = cell(R,1);
    for ii = 1:R
        ind = find(sourcemodelAAL.tissue == ii);
        voxc = mean(sourcemodel.pos(ind,:)); % centroid
        locs(ii,:) = voxc;        
        locsAAL{ii} = sourcemodel.pos(ind,:);       
    end
    
    if strcmp(roiopt,'g')
        locsc = locs;
        locs = cell2mat(locsAAL);
    end
end

%% Calculate lead fields

leadfield_name =sprintf( '%s/leadfields_%.0fmm.mat',processing_folder,gridres);
if ~exist(leadfield_name,'file')
    cfg                 = [];
    cfg.grad            = sens;
    cfg.headmodel       = vol;
    cfg.reducerank      = 2;
    cfg.channel         = {'MEG'};
    cfg.sourcemodel.pos = locs; %sourcemodel.pos
    cfg.sourcemodel.unit   = 'm';
    cfg.siunits         = true;
    cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
    [grid] = ft_prepare_leadfield(cfg);
    save(leadfield_name,'grid');
else
    load(leadfield_name);
end

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
icacomps = length(data.cfg.component);

C = cov(data_filt');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
%     Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*E(1)*eye(size(C)); % 1% max singular value =~ 70*noise

data.trial{1} = data_filt;
clear data_filt
%%
[~,~,ind]= intersect(sensall,data.label);
% whole evoked response time for sign uncertainty
cue_sample = cue_match.sample(cue_match.choice>0);
[datave,ttdel]= define_trials(cue_sample, data, tasktime, [.25,.4],0);
ntrials = length(datave.trial);
datavem = cell2mat(datave.trial);

trials =  cue_match.bv_index(cue_match.choice>0)-12;

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
gamble_cue = cue_match.choice(cue_match.choice>0);
gamble_cue(ttdel) = [];

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
sens_cue = squeeze(mean(datas(ind,:,:),2));


%% cue P300 
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

%% Choice

% whole evoked response time for sign uncertainty
choice_sample = choice_match.sample(choice_match.choice>0);
[datave,ttdel]= define_trials(choice_sample, data, tasktime, [.25,.4],0);
ntrials = length(datave.trial);
datavem = cell2mat(datave.trial);

trials =  choice_match.bv_index(choice_match.choice>0)-12;

% do mood and expectation (all 3 types) influence the P300?
% does not make too much sense to include RPE
E_choice = Es(trials);
E_choice(ttdel) = [];
Elta_choice = EltaH(trials);
Elta_choice(ttdel) = [];
Eltas_choice = E_LTA(trials);
Eltas_choice(ttdel) = [];

trials_choice = trials;

xi = choice_sample;
mood_choice = Fmood(xi); % interpolated mood timecourse
mood_choice(ttdel) = [];
Schoice= repmat(sub,length(mood_choice),1);
gamble_choice = choice_match.choice(choice_match.choice>0);
gamble_choice(ttdel) = [];

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
sens_choice = squeeze(mean(datas(ind,:,:),2));

%%
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

ltvchoice = table(Schoice,trials_choice',mood_choice',E_choice,Elta_choice,Eltas_choice,gamble_choice',...
    'VariableNames',{'subject','trial','mood','E','E_LTA','E_sum','choice'});

% save(save_name,'ltvcue','ltvchoice','sens_choice','sens_cue','Pcue');
save(save_name,'Pcue','-append');

end

return

%%
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/';
gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/mni_grid.txt');
ltv_cue = [];
ltv_choice = [];

S_cue = [];
S_choice = [];
P_cue = [];

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
    
    ltvchoice.recording = repmat(r,size(ltvchoice,1),1);
    ltv_choice = cat(1,ltv_choice,ltvchoice);
    
%     sens = ft_read_sens(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name],...
%         'senstype','meg');
%     [~,~,ind]= intersect(sensall,sens.label(strcmp(sens.chantype,'meggrad')));
            
    S_cue = cat(2,S_cue,sens_cue);   
    S_choice = cat(2,S_choice,sens_choice);
    P_cue = cat(2,P_cue,Pcue);
end

P_cue = P_cue(gridall==1,:);

dlmwrite([outpath,'Mcue_P300_30Hzlowpass.txt'],S_cue)
dlmwrite([outpath,'Mchoice_P300_30Hzlowpass.txt'],S_choice)
dlmwrite([outpath,'BFcue_P300_30Hzlowpass.txt'],P_cue)

writetable(ltv_cue,[outpath,'/latent_vars_cue.csv']);
writetable(ltv_choice,[outpath,'/latent_vars_choice.csv']);


%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';

if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat')
end

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/P300/';
freql = {'cue';'choice'};
for ff = 1:length(freql)
    latent_vars_name = sprintf('latent_vars_%s.csv',freql{ff});
    opts = detectImportOptions([data_path,latent_vars_name]);
    X = readtable([data_path,latent_vars_name],opts);
    fit_parameters = X.Properties.VariableNames(3:7);

    freq = sprintf('M%s_P300_30Hzlowpass',freql{ff});
    meg_data_name = sprintf('%s.txt',freq);
    meg = dlmread([data_path,meg_data_name]);

    outpath = sprintf('%s%s/',data_path,freq);
    nn =1;

    Tfit = cell(1,length(fit_parameters));
    pfit = cell(1,length(fit_parameters));
    for ii = 1:length(fit_parameters)
    cd([outpath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    Tfit{ii} = Xv.tStat;
    pfit{ii} = Xv.pValue;
    end

megz = zeros(size(meg));
for r = 1:18
    ind = X.recording == r;
    m = meg(:,ind);
    m = zscore(m(:));
    m = reshape(m,size(meg(:,ind)));
    megz(:,ind) = m;
end

T = struct;
T.label = sensall;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure(ff); clf; set(gcf,'color','w','position',[176 348 1262 385])

for ii = 1:length(fit_parameters)
subplot(1,5,ii)
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

figure(ff+2); clf
T.trial{1} = mean(megz,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-.4 .4];
ft_topoplotER(cfg, T)
saveas(gcf,sprintf('~/matlab/figures/%s_avg.png',freq))

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
saveas(gcf,sprintf('~/matlab/figures/%s_pvalue.png',freq))
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