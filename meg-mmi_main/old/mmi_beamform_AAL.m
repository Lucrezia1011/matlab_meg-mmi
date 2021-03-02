clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% AAL atlas 
gridres = 5; % resolution of beamformer grid in mm

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

%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

% Can do the same for ERF!
sn = 5;
    
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline

name_list = dir;
data_names = cell(1);
jj = 0;
for ii = 1:length(name_list)
    if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds')
        jj = jj+1;
        data_names{jj} = name_list(ii).name;
    end
end

% for runs = 1:length(data_names)
runs = 1;
% clearvars -except subn sn data_names sub data_path runs freq_band freq

data_name = data_names{runs};
%%


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

%% Sourcemodel warp MNI grid

% sourcemodel based on 5mm grid MNI brain
cfg = [];
cfg.mri = mri;
cfg.warpmni = 'yes';
cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
cfg.unit      = 'm';
cfg.nonlinear = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);

%% Find location of AAL ROIs
R = length(sourcemodelAAL.tissuelabel);
locs = zeros(R,3);
locsAAL = zeros(R,3);
for ii = 1:R
    ind = find(sourcemodelAAL.tissue == ii);
    locs(ii,:) = mean(sourcemodel.pos(ind,:));
    
    locsAAL(ii,:) = mean(template_grid.pos(ind,:));
end

figure;
ft_plot_headmodel(vol,'facecolor','none','edgecolor',[0.8 0.8 0.8]);
set(gcf,'color','w')
hold all
plot3(locs(:,1)*1000,locs(:,2)*1000,locs(:,3)*1000,'.')
%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.sourcemodel.pos = locs;
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
else
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = 'EEG';
    % cfg.demean = 'yes';
    % cfg.bpfilter = 'yes';
    % cfg.bpfreq = [1 150];
    try
        eog = ft_preprocessing(cfg);
        eog = eog.trial{1}(1,:);
    catch
        disp('Could not find EEG channel')
    end
    
    cfg =[];
    cfg.method = 'pca';
    comp_pca = ft_componentanalysis(cfg, data);
    score = comp_pca.trial{1}';
    compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
    icomp = nnz(compvar95) ;
    clc
    fprintf('%d components for 95perc. of data variance\n',icomp)
    
    if icomp>30
        disp('Reducing ICA components to 30')
        icomp = 30;
    end
    cfg =[];
    cfg.method = 'fastica';
    cfg.fastica.numOfIC = icomp;
    comp = ft_componentanalysis(cfg, data);
    
%     addpath('~/fieldtrip-20190812/forward/private/')

    figure
    cfg           = [];
    cfg.component = [1:icomp];       % specify the component(s) that should be plotted
    cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)
    
    
    cfg          = [];
    cfg.channel  = [1:5]; % components to be plotted
    cfg.viewmode = 'component';
    cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
    ft_databrowser(cfg, comp)
    
    try
        figure;
        plot(abs(corr(eog',comp.trial{1}')))
        xlabel('ICA component')
        ylabel('Correlation with EOG')
    end
    icadel = input('ICA component to eliminate (input as [''01'';''02'']): ');
    
    cfg = [];
    cfg.channel = cell(size(icadel,1),1);
    for ii = 1:size(icadel,1)
        cfg.channel{ii}  = ['fastica0',icadel(ii,:)];
    end
    
    [comps] = ft_selectdata(cfg, comp);
    save([processing_folder,'/ICA_artifacts'],'comps')
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data          = ft_rejectcomponent(cfg, comps,data);

%%

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[1 150],filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


%% Beamfomer
icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C); 
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

% Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value

L = grid.leadfield;

VE(1:size(L,2)) = {0};
W(1:size(L,2)) = {0};


for ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
   
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
%         if d(3) < 1
%             jj = 2; % The minumum singular value is degenerate
%         else
%             jj =3;
%         end
        lfo = lf*v(:,jj); % Lead field with selected orientation
%         w = Cr\lfo / (lfo'/Cr*lfo) ;  
        
%         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
%         Better Hall's or normalized weights
%         wnorm = w/sqrt(w'*noiseC*w);    
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;       
        W{ii} = w;
        VE{ii}  = w'*data.trial{1};
    
        clc
        fprintf('SAM running %d/%d .\n', ii, R)
      
end
clc
fprintf('Beamformer finished\n' )

VE = cell2mat(VE');

%%

VE_filt = ft_preproc_bandpassfilter(VE, data.fsample,[13 30],filt_order,'but');

newf = 150;
VE_filt = resample(VE_filt',newf,f);
% [ii,jj] = find(triu(ones(R),1)); % for multivariate leakge correction
[ii,jj] = find(~eye(R));

% VE_filt2 = zeros(size(VE_filt));
% VE_filt2(:,1:(108/2)) = VE_filt(:,1:2:107); 
% VE_filt2(:,(108/2 + 1):108) = VE_filt(:,2:2:108); 
% VE_filt2(:,109:end) = VE_filt(:,109:end);

AEC = zeros(R);

for k = 1:length(ii)
    x = leakage_reduction(VE_filt(:,ii(k)), VE_filt(:,jj(k)));
    ht = hilbert([x, VE_filt(:,jj(k))]) ;
    ht(1:30,:) = [];
    ht(end-30+1:end,:) = [];
    env = abs(ht);
    env = resample(env,1,newf); % downsample envelopes to 1Hz
    r = corr(env);
    
    AEC(ii(k),jj(k)) = r(1,2); 
    
    clc
    fprintf('FC running %.1f.\n', k/length(ii)*100)
    
%     figure; plot(VE_filt(:,ii(k)))
%     hold all
%     plot(VE_filt(:,jj(k)))
%     plot(x)
end

AEC = (AEC + AEC') /2;
figure; imagesc(AEC)

AECc = AEC(1:90,1:90);

% addpath('~/ppyll1/matlab/brain_map/ConGraph')
figure; go_3Dbrain_2019(AECc,locsAAL(1:90,:)*10,0.95);
axis equal
% 
% 
%% Read events
% 
bv_match = match_triggers_fc(data_name);

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;
buttonpress = bv_match.buttonpress;

%%

VE_filt = ft_preproc_bandpassfilter(VE, data.fsample,[4 8],filt_order,'but');

VE_fb = data;
VE_fb.label = sourcemodelAAL.tissuelabel;
VE_fb.trial{1} = VE_filt;

VE_fb = define_trials(outcome_match.sample(outcome_match.sample~=0),VE_fb,bv_match.time,[-3 4] );

cfg = [];
cfg.resamplefs = 300;
VE_fb = ft_resampledata(cfg,VE_fb);
IAC = zeros(R,R,length(VE_fb.trial{1}));

[ii,jj] = find(~eye(R));

ve = cell2mat(VE_fb.trial)';
for k = 1:length(ii)
    x = ve(:,ii(k));
    y = ve(:,jj(k));
    x = leakage_reduction(x, y);
    ht = hilbert([x, y]) ;
    
    % Eliminate edges from hilbert
    ht(1:30,:) = [];
    ht(end-30+1:end,:) = [];
    env = abs(ht);
    % calculate IAC on zscored envelopes
    envz = zscore(env);
    iac = envz(:,1).*envz(:,2);
    
    % Add back edges in order to reshape and average over trials
    iac(31:length(iac)+30) = iac;
    iac(end:end+30,:) = 0;
    iac = reshape(iac,length(VE_fb.trial{1}),size(VE_fb.trial,2));
        
    IAC(ii(k),jj(k),:) = mean(iac,2);
    
    clc
    fprintf('FC running %.1f.\n', k/length(ii)*100)
end

for k = 1:size(IAC,3)
    IAC(:,:,k) = (IAC(:,:,k) + IAC(:,:,k)' )/2;
end

time = linspace(-3,4,length(VE_fb.trial{1}));
for k = 0:.1:.4
    tl = k+[-0.1 0.1];
    twind = time>tl(1) & time < tl(2);
    figure
    go_3Dbrain_2019(mean(IAC(:,:,twind),3),locsAAL*10,0.98,[-.6 .6]);
end

figure; plot(time,squeeze(IAC(20,30,:)))
title(['IAC connectivity in \theta band: ',VE_fb.label{20},' and ',VE_fb.label{30}])
