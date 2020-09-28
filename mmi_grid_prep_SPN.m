% Stimulus Preceding negativity!

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

% Need this for evoked resposnse
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/';
% Y = dlmread([data_path,'/meg_trials_evoked_outcome.txt']);
% Y = reshape(Y, 360,116,size(Y,2));
% 
% opts = detectImportOptions([data_path,'latent_vars_evoked_outcome.csv']);
% Xv = readtable([data_path,'latent_vars_evoked_outcome.csv'],opts);

% function mmi_grid_prep_evoked(data_name,roiopt,gridres,freqband,filter_type)
% Anticipation period



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

% Need this for evoked responses
% Ym = Y(:,:,Xv.subject == str2double(sub));
% Ym = mean(Ym,3);
%% Read events

[bv_match,bv] = match_triggers_fc(data_name);

% cue_match = bv_match.answer;
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

% Standard model
A = bv.outcomeAmount; % Outcome
A(isnan(A)) = [];
ntrials = length(A);

Es = (bv.winAmount + bv.loseAmount )/2;
Es(isnan(Es)) = [];
RPEs = A - Es;

% LTA model
EltaH = cumsum(A)./(1:ntrials)'; % Expectation, defined by Hanna
EltaH(2:end+1) = EltaH;
EltaH(1) = 0; EltaH(end) = [];
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
EltaH(2:end+1) = EltaH;
EltaH(1) = (bv.winAmount(13)+bv.loseAmount(13))/2;

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
    %% Find location of AAL ROIs
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

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,8,filt_order,'but');
icacomps = length(data.cfg.component);

C = cov(data_filt');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
%     Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*E(1)*eye(size(C)); % 1% max singular value =~ 70*noise

data.trial{1} = data_filt;
clear data_filt

dataf = data;
dataf.trial{1} = data_filt;
gamble_sample = choice_match.sample(choice_match.gamble==1);
% whole evoked response time for sign uncertainty
[datave,ttdel]= define_trials(gamble_sample, dataf, tasktime, [-1,5],0);
trials =  choice_match.bv_index(choice_match.gamble==1)-12;
trials(ttdel) = [];
rpe =RPEs(trials); 
E = Es(trials);
certain_sample = choice_match.sample(choice_match.choice==3);
% whole evoked response time for sign uncertainty
[datave1,ttdel]= define_trials(certain_sample, dataf, tasktime, [-1,5],0);

timew = 3:.2:3.8; clf
for c = 1:4
for ii = 1:length(timew)
subplot(4,length(timew),ii + (c-1)*length(timew))
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-2 2]*1e-14;
cfg.xlim = timew(ii) + [0 0.2];
cfg.baseline = [2 3];
switch c
    case 1
%         cfg.trials = find(trials<=27); 
        cfg.trials = find(abs(rpe)>5); 
        ft_topoplotER(cfg, datave)
    case 2
%         cfg.trials = find(trials>27 & trials<=54); 
        cfg.trials = find(abs(rpe)<5); 
        ft_topoplotER(cfg, datave)
    case 3
%         cfg.trials = find(trials>54); 
        cfg.trials = find(abs(rpe)>5 & E<2); 
        ft_topoplotER(cfg, datave)
    case 4
        cfg.trials = 'all'; 
        ft_topoplotER(cfg, datave1)
end
drawnow
% title('RPE LTA')
end
end


ntrials = length(datave.trial);
datavem = cell2mat(datave.trial);

trials =  choice_match.bv_index(choice_match.gamble==1)-12;
trials(ttdel) = [];

expectation = Es(trials);
% rpe = RltaH(trials);

EltaH = EltaH(trials);

E_LTA = E_LTA(trials);

xi = gamble_sample;
mood = Fmood(xi); % interpolated mood timecourse
mood(ttdel) = [];
S = repmat(sub,length(mood),1);

%% Sensor level

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);

Sant = squeeze(mean(abs(datas),2));

% timew = datave.time{1}>=.100 & datave.time{1}<=.150;  % ACC in 23911
% S100150 =  squeeze(mean(datas(:,timew,:),2));
% 
% timew = datave.time{1}>=.175 & datave.time{1}<=.225;  % ACC in 23911
% S175225 =  squeeze(mean(datas(:,timew,:),2)); %P200
% 
% timew = datave.time{1}>=.225 & datave.time{1}<=.275;  % ACC in 23911
% S225275 =  squeeze(mean(datas(:,timew,:),2));  % largest peak
% 
% timew = datave.time{1}>=.250 & datave.time{1}<=.350;  % ACC in 23911
% S250350 =  squeeze(mean(datas(:,timew,:),2));  % FRN for Hajack
% 
% timew = datave.time{1}>=.300 & datave.time{1}<=.500;  % ACC in 23911
% S300500 =  squeeze(mean(datas(:,timew,:),2));  % P300


%% Beamfomer

L = grid.leadfield(grid.inside);

Pant = cell(size(L));
% P100150 = cell(size(L));
% P175225 = cell(size(L));
% P225275 = cell(size(L));
% P250350 = cell(size(L));
% P300500 = cell(size(L));
for ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    
    % %  G O'Neill method, equivalent to ft
    [v,d] = svd(lf'/Cr*lf);
    d = diag(d);
    jj = 2;
    
    lfo = lf*v(:,jj); % Lead field with selected orientation
    
    w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
    
    ve = w'*datavem;
    ve = reshape(ve,[datave.sampleinfo(1,2),ntrials]);  
   
    Pant{ii} = mean(abs(ve),1);
    
%     timew = datave.time{1}>=.100 & datave.time{1}<=.150;  % ACC in 23911
%     P100150{ii} = mean(ve(timew,:),1)*sign(c(ind));
%     
%     timew = datave.time{1}>=.175 & datave.time{1}<=.225;  % ACC in 23911
%     P175225{ii} = mean(ve(timew,:),1)*sign(c(ind)); %P200
%     
%     timew = datave.time{1}>=.225 & datave.time{1}<=.275;  % ACC in 23911
%     P225275{ii} = mean(ve(timew,:),1)*sign(c(ind));  % largest peak
%     
%     timew = datave.time{1}>=.250 & datave.time{1}<=.350;  % ACC in 23911
%     P250350{ii} =  mean(ve(timew,:),1)*sign(c(ind));  % FRN for Hajack
%     
%     timew = datave.time{1}>=.300 & datave.time{1}<=.500;  % ACC in 23911
%     P300500{ii} =   mean(ve(timew,:),1)*sign(c(ind));  % P300
%     
    if mod(ii,100) == 0
        clc
        fprintf('%s\nSAM running %.1f\n',...
            data_name,ii/length(L)*100)
    end

end

Pin = cell2mat(Pant');
Pant = zeros(size(grid.leadfield,2),size(Pin,2));
Pant(grid.inside,:) = Pin;

% Pin = cell2mat(P175225');
% P175225 = zeros(size(grid.leadfield,2),size(Pin,2));
% P175225(grid.inside,:) = Pin;
% 
% Pin = cell2mat(P200300');
% P200300 = zeros(size(grid.leadfield,2),size(Pin,2));
% P200300(grid.inside,:) = Pin;
% 
% Pin = cell2mat(P225275');
% P225275 = zeros(size(grid.leadfield,2),size(Pin,2));
% P225275(grid.inside,:) = Pin;
% 
% Pin = cell2mat(P250350');
% P250350 = zeros(size(grid.leadfield,2),size(Pin,2));
% P250350(grid.inside,:) = Pin;
% 
% Pin = cell2mat(P300500');
% P300500 = zeros(size(grid.leadfield,2),size(Pin,2));
% P300500(grid.inside,:) = Pin;

%%
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/';
save_name = sprintf('%s%s',outpath,sub);

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end


ltvant = table(S,trials',mood',expectation,E_LTA,EltaH,'VariableNames',...
    {'subject','trial','mood','E','E_sum','E_LTA'});

%%
save(save_name,'ltvant','Pant','Sant','-append');

end

return
%% Find all common gradiometers
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/';

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
    sens = ft_read_sens(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name],...
        'senstype','meg');
    if s == 1
        sensall = sens.label(strcmp(sens.chantype,'meggrad'));
    else
        sensall = intersect(sensall,sens.label(strcmp(sens.chantype,'meggrad')));
    end
end
%%
ltv = [];
P= [];
% p175225= [];
% p200300= [];
% p250350= [];
% p300500= [];
% p225275= [];

S= [];
% s175225= [];
% s200300= [];
% s250350= [];
% s300500= [];
% s225275= [];

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
    ltvant.recording = repmat(r,size(ltvant,1),1);
    ltv = cat(1,ltv,ltvant);
    
    sens = ft_read_sens(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name],...
        'senstype','meg');
    [~,~,ind]= intersect(sensall,sens.label(strcmp(sens.chantype,'meggrad')));
   
         
    S = cat(2,S,Sant(ind,:));
%     s175225 = cat(2,s175225,S175225(ind,:));
%     s200300 = cat(2,s200300,S200300(ind,:));
%     s225275 = cat(2,s225275,S225275(ind,:));
%     s250350 = cat(2,s250350,S250350(ind,:));
%     s300500 = cat(2,s300500,S300500(ind,:));
   
    P = cat(2,P,Pant);
%     p175225 = cat(2,p175225,P175225);
%     p200300 = cat(2,p200300,P200300);
%     p225275 = cat(2,p225275,P225275);
%     p250350 = cat(2,p250350,P250350);
%     p300500 = cat(2,p300500,P300500);
end

gridall = dlmread([outpath,'mni_grid.txt']);
P = P(gridall==1,:);
% p175225 = p175225(gridall==1,:);
% p200300 = p200300(gridall==1,:);
% p250350 = p250350(gridall==1,:);
% p300500 = p300500(gridall==1,:);
% p225275 = p225275(gridall==1,:);


dlmwrite([outpath,'BF_anticipation_30Hzlowpass.txt'],P)
% dlmwrite([outpath,'BF_100-150ms_30Hzlowpass.txt'],p100150)
% dlmwrite([outpath,'BF_175-225ms_30Hzlowpass.txt'],p175225)
% dlmwrite([outpath,'BF_225-275ms_30Hzlowpass.txt'],p225275)
% dlmwrite([outpath,'BF_250-350ms_30Hzlowpass.txt'],p250350)
% dlmwrite([outpath,'BF_300-500ms_30Hzlowpass.txt'],p300500)

dlmwrite([outpath,'M_anticipation_30Hzlowpass.txt'],S)
% dlmwrite([outpath,'M_100-150ms_30Hzlowpass.txt'],s100150)
% dlmwrite([outpath,'M_175-225ms_30Hzlowpass.txt'],s175225)
% dlmwrite([outpath,'M_225-275ms_30Hzlowpass.txt'],s225275)
% dlmwrite([outpath,'M_250-350ms_30Hzlowpass.txt'],s250350)
% dlmwrite([outpath,'M_300-500ms_30Hzlowpass.txt'],s300500)

writetable(ltv,[outpath,'/latent_vars_anticipation.csv']);

%%


Sm = mean(S,1); 
clc
for ii = 3:6
fit_parameter = ltv.Properties.VariableNames{ii};
lme_formula = sprintf('MEG ~ %s +(-1 + %s|subject) + (1|trial) + (1|recording)',fit_parameter,fit_parameter);
ltv.MEG = Sm';
lme = fitlme(ltv,lme_formula); %Fixed effects for RPE
fprintf('%s: p-value = %.4f\n',fit_parameter,lme.Coefficients.pValue(2));
end

%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';
times = {'100-150';'175-225';'200-300';'225-275';'250-350';'300-500'};
for ii = 1:6
datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/M_%sms_30Hzlowpass/',times{ii});
nn =1;
cd([datapath,'lme_RPE'])
opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
Trpe{1} = Xv.tStat;

cd([datapath,'lme_RPE_abs'])
opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
Trisk{1} = Xv.tStat;

cd([datapath,'lme_RPE_LTA'])
opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
Tlta{1} = Xv.tStat;

T = struct;
T.label = sensall;
T.time{1} = str2double(times{ii}(1:3));
T.trial = Tlta;
T.sampleinfo = [1 1];

figure(ii); clf; set(gcf,'color','w','position',[176 348 1253 494])
subplot(131)
T.avg = Trpe; T.trial = Trpe;
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-4 4];
ft_topoplotER(cfg, T)
title('RPE')

subplot(132)
T.avg = Trisk; T.trial = Trisk;
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-4 4];
ft_topoplotER(cfg, T)
title('|RPE|')

subplot(133)
T.avg = Tlta; T.trial = Tlta;
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-4 4];
ft_topoplotER(cfg, T)
title('RPE LTA')

saveas(gcf,sprintf('~/matlab/figures/M_%sms_30Hzlowpass.png',times{ii}))
end
%% Plot Linear mixed effects model for grid
times = {'100-150';'175-225';'200-300';'225-275';'250-350';'300-500'};

ii = 4%[3,4,6];
mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/mni_grid.txt');
Tlta(1) = {zeros(size(gridall))};
Trisk(1) = {zeros(size(gridall))};
Trpe(1) = {zeros(size(gridall))};

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


datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/BF_%sms_30Hzlowpass/',times{ii});
    cd([datapath,'lme_RPE_LTA'])
    X = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    end
    Tlta{1}(gridall==1) = X;
    
    
    cd([datapath,'lme_RPE'])
    
    X = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    end
    
    Trpe{1}(gridall==1) = X;
    
    
    cd([datapath,'lme_RPE_abs'])
    
    X = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    end
    Trisk{1}(gridall==1) = X;
    
    
    fprintf(sprintf('Loaded'))
    

%%
% close all
freq = 1;
% Te{freq}(Te{freq}==0) = 1;
% sourceant.pow = log10(Te{freq});
for p = 1:3
    sourceant =[];

    switch p
        case 1
            sourceant.pow = Trpe{freq};
            fit_parameter = 'RPE';
        case 2
            sourceant.pow = Trisk{freq};
            fit_parameter = '|RPE|';
        case 3
            sourceant.pow = Tlta{freq};
            fit_parameter = 'RPE LTA';
    end
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
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
title(sprintf('%s \npeak t-value %.1f',...
    fit_parameter,max(abs(sourceant.pow(:)))))
if strcmp(cfg.method,'slice')
saveas(gcf,sprintf('~/matlab/figures/BF_%sms_30Hzlowpass_%.f.0.png',times{ii},p))
end
end