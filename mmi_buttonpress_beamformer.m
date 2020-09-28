function sourcenorm = mmi_buttonpress_beamformer(data_name,freqband,wl,mu)

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% Co-register MRI from fiducial positions

% LTA model latent variables:
% EC: Expectation of certain value
% EG: Expectation during gabling
% Ediff: Drift rate
% LTA: Long term average with gamma:   1/t * sum_i=1 ^t(V(i)^gamma),   cumsum(LTA.OutcomeAmount^gamma)./(1:ntrials)'
% V_i^gamma = outcome of trial i
% new_p = subjective winning probability
% RPE = Reward prediction error
% LTA_sum  = sum(LTA)
% RPE_sum = sum(RPE)
% log_like
% mood_log_like

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)


bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
for ii = 1:length(bv_names)
    if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
        bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
ltv = readtable(bv_name,opts); % bahavioral data
%% Co-register MRI

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

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
gridres = 5; % resolution of beamformer grid in mm

% Load fieldtrip 10mm MNI grid
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
template_grid = sourcemodel;

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

%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.sourcemodel.pos = sourcemodel.pos;
cfg.sourcemodel.unit   = 'm';
cfg.siunits         = true;
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);

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
%     data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[1 150],filt_order,'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,50,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt

%% Read events

bv_match = match_triggers_fc(data_name);

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;

%% Filter data in theta band
% freqband = [4,8];

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, freqband,filt_order,'but')';
% ICA clean reduces min(svd(C)) by several orders of magnitude while
% leaving 
C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));

% Cr = C + mu*noise;
Cr = C + mu*max(svd(C)).*eye(size(C));


%% Button presses 

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = {'UADC005';'UADC006';'UADC007'}; % Time channel!
buttons = ft_preprocessing(cfg);

buttonsd = diff(buttons.trial{1}');
buttonpress = buttonsd>1.5;
[samples,~] = find(buttonpress);
samples = sort(samples);

samples(find(diff(samples)<0.2*f)+1) = [];


wind_press = 1:wl*f; 

% cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
if samples(end)+wind_press(end) > data.sampleinfo(2)
    samples(end) = [];
end

Cah = zeros([size(C),length(samples)]);
for tt= 1:length(samples)
    data_win = data_filt(samples(tt)+wind_press,:);
    Cah(:,:,tt) = cov(data_win);
end

ind = ITI_match.sample>0;
samples_iti = ITI_match.sample(ind);
wind_press = 1:2.0*f;

if samples_iti(end)+wind_press(end) > data.sampleinfo(2)
    samples_iti(end) = [];
end

r = floor(2/wl);

Cal = zeros([size(C),size(samples_iti,2)*r]);
for tt= 1:size(samples_iti,2)
    data_lose = data_filt(samples_iti(tt)+wind_press,:);
%     Cal(:,:,tt) = cov(data_lose);
    for tti = 1:r
        Cal(:,:,(tt-1)*r+tti) = cov(data_lose((wl*f*(tti-1))+1:wl*f*tti,:));
    end
end

if length(samples)>(r*length(samples_iti))
    Cah = Cah(:,:,1:(r*length(samples_iti)));
elseif length(samples)<(r*length(samples_iti))
    Cal = Cal(:,:,randperm(r*length(samples_iti),length(samples)));
%     Cal = Cal(:,:,1:length(samples));
end
Cah = mean(Cah,3);
Cal = mean(Cal,3);


%% Mood samples
% ind = mood_match.sample>0;
% mood_match.mood(~ind)= NaN;
% 
% v1 = mood_match.mood(ind);
% samples1 = mood_match.sample(ind);
% 
% ind = blockmood_match.sample>0;
% blockmood_match.mood(~ind) = NaN;
% v0 = blockmood_match.mood(ind);
% samples0 = blockmood_match.sample(ind);
% 
% v = [v1,v0];
% [s,ind] = sort(v);
% 
% samples = [samples1, samples0];
% samples = samples(ind);
% 
% hs = s > median(s);
% ls = s < median(s);
% 
% happy_samples = samples(hs);
% sad_samples = samples(ls);
% 
% if nnz(hs)> nnz(ls)
%     happy_samples(1) = [];
% elseif nnz(hs)< nnz(ls)
%     sad_samples(end) = [];
% end
% 
% wind = 0*f:3*f; % Exclude edges for sensory and motor overlap
% % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
% 
% 
% Cah = zeros([size(C),size(happy_samples,2)]);
% for tt= 1:size(happy_samples,2)
%     data_win = data_filt(happy_samples(tt)+wind,:);
%     Cah(:,:,ii) = cov(data_win);
% end
% Cah = mean(Cah,3);
% 
% Cal = zeros([size(C),size(sad_samples,2)]);
% for tt= 1:size(sad_samples,2)
%     data_lose = data_filt(sad_samples(tt)+wind,:);
%     Cal(:,:,ii) = cov(data_lose);
% end
% Cal = mean(Cal,3);


%% Beamformer
L = grid.leadfield;
W = cell(size(L));

Zstat_mood(1:size(L,2)) = {0};
Tstat_mood(1:size(L,2)) = {0};

for ii = 1:length(L)
    lf = L{ii};
    if ~isempty(lf)
        
        % %           G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        if d(3) < 1
            jj = 2; % The minumum singular value is degenerate
        else
            jj =3;
        end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        w = Cr\lfo / (lfo'/Cr*lfo) ;
        
        %             Qap = w'*Cap*w;
        %             Qan = w'*Can*w;
        
        Qah = w'*Cah*w;
        Qal = w'*Cal*w;
        
        n = w'*noise*w;
        
        W{ii} = w;
        %             Tstat{ii} = (Qa - Qc) ./ (na + nc);
        %             Tstat_rpe{ii} = (Qap - Qan) ./ (Qap + Qan); % normalised version
        Tstat_mood{ii} = (Qah - Qal) ./ (Qah + Qal); % normalised version
%         Zstat_mood{ii} = ((Qah + Qal)/2 ) ./ n;
        
    end
    
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
    
end
clc
fprintf('SAM finished\n' )
%% Save result

if ~exist([data_path,'results'],'dir')
    cd(data_path)
    mkdir results
end
cd([data_path,'results'])

mriname = 'coreg_anat';
if ~exist('coreg_anat_ortho+orig.BRIK','file')
    writebrik(mriname,mri)
end

cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';

%     if ~exist(zname,'file')

% sourceTstat.avg.pow =  cell2mat(Zstat_mood);
% sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
% sourcePostInt.anatomy = sourcePostInt.pow;

% zname =  sprintf('MoodZ_3s_PseudoT_%d-%dHz_mu%.0s',freqband(1),freqband(2),mu);
% zname_nii = [zname,'.nii'];
% ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');


sourceTstat.avg.pow =  cell2mat(Tstat_mood);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

% 3s before mood rating high - low mood
zname =  sprintf('Buttonpress_%.1fs_PseudoT_%d-%dHz_mu%.0s',wl,freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
writebrik(zname,sourcePostInt)
% ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');

%%
%     crang = [];
%     sourcePostInt.anatomy = mri.anatomy;
%     cfg = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter = 'pow';
%     cfg.maskparameter = 'pow';
%     cfg.funcolormap  = 'auto';
%     cfg.funcolorlim   = crang;
%     cfg.opacitylim = crang;
%     ft_sourceplot(cfg, sourcePostInt);

%% Template transform

sourcePostInt.anatomy = mri.anatomy;
cfg = [];
cfg.parameter = 'pow';
cfg.template = '/home/liuzzil2/MNI152_T1_2009c.nii';
[sourcenorm] = ft_volumenormalise(cfg, sourcePostInt);

% if ~exist('ft_coreg_anat_norm.nii','file')
%     zname =  sprintf('ft_coreg_anat_norm');
%     zname_nii = [zname,'.nii'];
%     ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');
% end

sourcenorm.anatomy = sourcenorm.pow;
zname =  sprintf('Buttonpress_%.1fs_PseudoT_%d-%dHz_norm_mu%.0s',wl,freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');


% sourceTstat = ft_read_mri(sprintf('Choicepress_0.5s_PseudoT_%d-%dHz_mu%.0s.nii',freqband(1),freqband(2),mu));
% sourceTstat.pow = sourceTstat.anatomy;
% sourceTstat.pow(isnan(sourceTstat.pow)) = 0;
% sourceTstat.anatomy = mri.anatomy;
% sourceTstat.brain = segmentmri.brain;
% cfg = [];
% cfg.method = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.maskparameter = 'pow';
% 
% cfg.funcolorlim = 'maxabs';
% cfg.funcolorlim = [-.5 .5];
% ft_sourceplot(cfg,sourceTstat)

%%

% sourcePostInt.pow(isnan(sourcePostInt.pow)) = 0;
% 
% cfg = [];
% cfg.method = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.maskparameter = 'pow';
% 
% cfg.funcolorlim = 'maxabs';
% cfg.funcolorlim = [];
% ft_sourceplot(cfg,sourcePostInt)



