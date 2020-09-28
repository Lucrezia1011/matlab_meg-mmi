function mmi_grid_prep_PowerITI(data_name,roiopt,gridres,freqband)
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
% roiopt = 's' sensors
% roiopt = 'grid' mni grid
% gridres = grid resolution in mm, for 'g' and 'grid' options

% addpath /home/liuzzil2/fieldtrip-20190812/
% ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca

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

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

%% Read events

[bv_match,bv] = match_triggers_fc(data_name);

% cue_match = bv_match.answer;
% choice_match = bv_match.choice;
% outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;

ITI_sample = bv_match.ITI.sample(bv_match.ITI.sample~=0); 

% ITI_sample0 = bv_match.ITI.sample(bv_match.ITI.bv_index~=bv_match.ratemood.bv_index); 
% ITI_sampler = bv_match.ITI.sample(bv_match.ITI.bv_index==bv_match.ratemood.bv_index & bv_match.ITI.sample~=0); 

indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;
[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

mood = Fmood(ITI_sample); % interpolated mood timecourse

trials =  bv_match.ITI.bv_index(bv_match.ITI.sample~=0) -12;
% trials = cat(2,trials,bv_match.blockmood.bv_index(bv_match.blockmood.sample~=0)-0.5);

% Standard model
A = bv.outcomeAmount; % Outcome
A(isnan(A)) = [];
ntrials = length(A);
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

E_LTA = E_LTA(trials);
RPE = RPE(trials);

Elta = EltaH(trials);
Rlta = RltaH(trials);
%%

if strcmp(roiopt,'sens')
    
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
    
    filt_order = []; % default
    
    data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,freqband,filt_order,'but');
    
    data.trial{1} = data_filt;
    clear data_filt
    
    [datave,ttdel]= define_trials(ITI_sample, data, tasktime, [0,2],0);
    ntrials = length(datave.trial);
    %% Sensor level
    datavem = cell2mat(datave.trial);
    datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
    
    V = squeeze(var(datas,0,2));
    
    
    save_name = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/ITI/%.0f-%.0fHz_%s',...
        freqband(1),freqband(2),sub);

    n = str2double(data_name(end-3));
    if ~isnan(n) %check for number at end of filename
        save_name = [save_name,'_',data_name(end-3)];
    else
        save_name = [save_name,'_1'];
    end

    mood(ttdel) = [];
    trials(ttdel) = [];
    S = repmat(sub,length(mood),1);
    RPE(ttdel) = [];
    E_LTA(ttdel) = [];
    Rlta(ttdel) = [];
    Elta(ttdel)  = [];
    ltvmood = table(S,trials',mood',Elta,E_LTA,Rlta,RPE,'VariableNames',...
        {'subject','trial','mood','E_LTA','E_sum','RPE_LTA','RPE_sum'});
    
    save(save_name,'V','ltvmood');

else
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

data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,freqband,filt_order,'but');
icacomps = length(data.cfg.component);

C = cov(data_filt');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
%     Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*E(1)*eye(size(C)); % 1% max singular value =~ 70*noise

data.trial{1} = data_filt;
clear data_filt

[datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
ntrials = length(datave.trial);
Ctt = zeros([size(C),ntrials]);
for tt=1:ntrials
    Ctt(:,:,tt) = cov(datave.trial{tt}');
end


%% Beamfomer

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
    
    pp  = zeros(ntrials,1);
    for tt = 1:ntrials
        pp(tt) =  w'*Ctt(:,:,tt)*w;
    end
    P{ii} = pp;
    
    if mod(ii,100) == 0
        clc
        fprintf('%s\n%.0f-%.0fHz: SAM running %.1f\n',...
            data_name,freqband(1),freqband(2),ii/length(L)*100)
    end

end

P  = cell2mat(P)';

%%

save_name = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/%.0f-%.0fHz_%s',...
    freqband(1),freqband(2),sub);

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end

mood(ttdel) = [];
trials(ttdel) = [];
S = repmat(sub,length(mood),1);
RPE(ttdel) = [];
E_LTA(ttdel) = [];
ltvmood = table(S,trials',mood',RPE,E_LTA,'VariableNames',...
    {'subject','trial','mood','RPE','E'});

%%
save(save_name,'ltvmood','P');

end






