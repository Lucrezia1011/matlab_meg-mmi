function grid = mniLeadfields(data_name,processing_folder,gridres,mri)
% Calculate lead fields with fieldtrip from Nolte's single shell    
% approximation on regularly spaced MNI grid warped to individual anatomy 
% 
% grid = mniLeadfields(data_name,processing_folder,gridres,mri)
% data_name          = dataset name (.ds)
% processing_folder  = folder for data derivatives
% gridres            = beamformer grid resolution in mm
% mri                = co-registered mri

leadfield_name =sprintf( '%s/leadfields_%.0fmm.mat',processing_folder,gridres);
if ~exist(leadfield_name,'file')
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
    
   
    %% MNI template brain
    
    % Load fieldtrip MNI grid
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
    locs = sourcemodel.pos;
   
    %% Calculate lead fields
      
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
    
    %% Eliminate Bad channels    
   
    % Get Bad channel names
    fid = fopen([data_name,'/BadChannels']);
    BadChannels = textscan(fid,'%s');
    fclose(fid);

    % Delete Bad channels
    chanInd = zeros(size(grid.label));
    for iiC = 1:length(BadChannels{1})
        chanInd = chanInd | strcmp(grid.label,BadChannels{1}{iiC});
    end
    grid.label(find(chanInd)) = [];
    for ii = find(grid.inside)'
        grid.leadfield{ii}((find(chanInd)),:) = [];
    end
    %%
    
    save(leadfield_name,'grid');
else
    load(leadfield_name);
end
