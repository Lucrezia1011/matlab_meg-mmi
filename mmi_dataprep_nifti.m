function dataprep = mmi_dataprep_nifti(sub)

clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

sub = '24138';
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline

if ~exist(data_name,'dir')
    data_name = [sub,'MMI_mmi3_proc1.ds'];
end

if strcmp(sub, '24201') || strcmp(sub, '22695')
    data_name = [sub,'MMI_mmi3_proc2.ds'];
end

% if exist([data_path,'results/ft_coreg_anat.nii'],'file')
%     mri = ft_read_mri([data_path,'results/ft_coreg_anat.nii']);
%     mri.coordsys = 'ctf';
% else
    mri_name = [sub,'_anat.nii'];
    
    if ~exist(mri_name,'file')
        unix(['gunzip ',mri_name])
    end
    
    mri = ft_read_mri(mri_name,'dataformat','nifti');
    fidtag = fopen('fids.tag');
    [firsthead, pos] = textscan(fidtag,"%s", 7); 
    fids = textscan(fidtag,'%q%f%f%f%u%u');
    fclose(fidtag);
    
    tagset_coord = zeros(3,3);
    tagset_coord(:,1) = fids{2};
    tagset_coord(:,2) = fids{3};
    tagset_coord(:,3) = fids{4};
    
    
    m = [   -1  0   0  mri.dim(1)
        0   -1   0   201
        0   0   1   0
        0   0   0   1] ;
    
%     mri.transform(3,3) = 1;
%     fiducial_coordn = (mrin.hdr.vox2ras0 \[tagset_coordn,ones(3,1)]')';
%     fiducial_coordn = [tagset_coordn,ones(3,1)]/(mrin.transform )';
    fiducial_coord = ((mri.hdr.vox2ras0/m) \[tagset_coord,ones(3,1)]')';
    
    %%
    
    figure(1); clf;
%     subplot(131)
    imagesc(squeeze(mri.anatomy(154,:,:))')
    set(gca,'Ydir','default','Xdir','reverse'); axis equal; colormap gray
    figure(2); clf
%     subplot(132)
    imagesc(squeeze(mri.anatomy(:,139,:))')
    set(gca,'Ydir','default'); axis equal;colormap gray; caxis([0 700])
%     subplot(133)
figure(3); clf
    imagesc(squeeze(mri.anatomy(:,:,88))')
    set(gca,'Ydir','default'); axis equal
    colormap gray
    %%
    figure; subplot(231)
    imagesc(squeeze(mri.anatomy(153,:,:)))
    subplot(232)
    imagesc(squeeze(mri.anatomy(:,139,:)))
    subplot(233)
    imagesc(squeeze(mri.anatomy(:,:,88)))
    
    subplot(234)
    imagesc(squeeze(mri.anatomy(100,:,:)))
    subplot(235)
    imagesc(squeeze(mri.anatomy(:,150,:)))
    subplot(236)
    imagesc(squeeze(mri.anatomy(:,:,150)))
    %%
    mri_name = [sub,'_anat+orig.BRIK'];
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
    cfg.viewresult = 'yes';
    
    mri = ft_volumerealign(cfg,mri);
    
%     ft_write_mri([data_path,'results/ft_coreg_anat.nii'],mri,'dataformat','nifti');
% end

%% Segment MRI
cfg = [];
cfg.output  = 'brain';
segmentmri = ft_volumesegment(cfg,mri);

%% Head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentmri);
sens = ft_read_sens(data_name,'senstype','meg');

%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.resolution      = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit   = 'cm';
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)

if ~exist([data_path,'/lead_fields.mat'],'file')
    [grid] = ft_prepare_leadfield(cfg);
    save([data_path,'/lead_fields'],'grid')
else
    load([data_path,'/lead_fields.mat'])
end
%% Clean Data

% Unsure of best pre-processing for baseline in sub 24071
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
% cfg.demean = 'yes';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [1 150];
data = ft_preprocessing(cfg);
f = data.fsample;

if exist([data_path,'ICA_artifacts.mat'],'file')
    load([data_path,'ICA_artifacts.mat']);
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
    save([data_path,'/ICA_artifacts'],'comps')
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data_clean    = ft_rejectcomponent(cfg, comps,data);

%% Read events

bv_match = match_triggers_fc(data_name);

%%
dataprep = [];
dataprep.data = data_clean;
dataprep.mri = mri;
dataprep.leadfield = grid;
dataprep.bv = bv_match;
dataprep.dataname = data_name;

