function dataprep = mmi_dataprep_mni(data_name,gridres,mnigrid,icaopt)
% Specify subject SDAN and option for warped mni grid (yes,no). 
% gridres = beamformer grid resolution in mm

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

% if exist([data_path,'results/ft_coreg_anat.nii'],'file')
%     mri = ft_read_mri([data_path,'results/ft_coreg_anat.nii']);
%     mri.coordsys = 'ctf';
% else
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
    
%     ft_write_mri([data_path,'results/ft_coreg_anat.nii'],mri,'dataformat','nifti');
% end

%% Segment MRI
if ~exist([processing_folder,'/headmodel.mat'],'file')
    cfg = [];
    cfg.output  = 'brain';
    segmentmri = ft_volumesegment(cfg,mri);

%% Head model

    cfg = [];
    cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, segmentmri);
    
    save([processing_folder,'/headmodel.mat'],'vol')
else
    load([processing_folder,'/headmodel.mat']);
end
sens = ft_read_sens(data_name,'senstype','meg');

%% Sourcemodel warp MNI grid

if mnigrid == true
    % Load fieldtrip 10mm MNI grid
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    template_grid = sourcemodel;
    clear sourcemodel
    
    cfg = [];
    cfg.mri = mri;
    cfg.warpmni = 'yes';
    cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
    cfg.unit      = 'm';
    cfg.nonlinear = 'yes';
    sourcemodel = ft_prepare_sourcemodel(cfg);
      
else
    
    cfg = [];
    cfg.headmodel  = vol;
    cfg.unit     = 'm';
    cfg.resolution  = gridres/1000; % convert to meters
    sourcemodel = ft_prepare_sourcemodel(cfg);
end
% visualise normalized grid on individual brain model
% figure; hold on; set(gcf,'color','w')
% ft_plot_headmodel(vol, 'edgecolor','none','facealpha',0.4,'facecolor',[0.5 0.5 0.5])
% ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)*1000) % *1000 as vol is expressed in mm


% mnitemplate = ft_read_mri('/home/liuzzil2/MNI152_T1_2009c.nii','dataformat','nifti');


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

% if ~exist([data_path,'/lead_fields.mat'],'file')
%     [grid] = ft_prepare_leadfield(cfg);
%     save([data_path,'/lead_fields'],'grid')
% else
%     load([data_path,'/lead_fields.mat'])
% end

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

if icaopt == true

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
data_clean    = ft_rejectcomponent(cfg, comps,data);
else
    data_clean = data;
    data_clean.cfg.component = 0;
end
% 
% %% Read events
% 
% bv_match = match_triggers_fc(data_name);

%%
dataprep = [];
dataprep.data = data_clean;
dataprep.mri = mri;
dataprep.leadfield = grid;
dataprep.headmodel = vol;
dataprep.sourcemodel = sourcemodel;
% dataprep.bv = bv_match;
dataprep.dataname = data_name;
