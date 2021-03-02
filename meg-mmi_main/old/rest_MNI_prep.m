clear all
close all
clc

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn =1:16 %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
    
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_rest_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = cell(1);
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18)
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    if jj>0
        for runs = 1:length(data_names)
            zz = zz +1;
            param_list{zz} = data_names{runs};
        end
    end
end

% function rest_MNI_prep(data_name)
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%%
for zz = [1:length(param_list)]
    %% Co-register MRI from fiducial positions
    data_name = param_list{zz};
    sub = data_name(1:5);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
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
    
    
    %
    % if strcmp(roiopt,'g'
    %     locsc = locs;
    %     locs = cell2mat(locsAAL);
    % end
    %% Calculate lead fields
    if ~exist([processing_folder,'/leadfields_5mm_MNI.mat'],'file')
        
        cfg                 = [];
        cfg.grad            = sens;
        cfg.headmodel       = vol;
        cfg.reducerank      = 2;
        cfg.channel         = {'MEG'};
        % cfg.sourcemodel.pos = locs;
        cfg.sourcemodel.pos = sourcemodel.pos;
        cfg.sourcemodel.unit   = 'm';
        cfg.siunits         = true;
        cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
        grid = ft_prepare_leadfield(cfg);
        save([processing_folder,'/leadfields_5mm_MNI.mat'],'grid');
    else
        load([processing_folder,'/leadfields_5mm_MNI.mat']);
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
    
    
    % need to automatically detect electrical artefacts, high variance noise
    % and muscle noise (for beta and gamma).
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
        
        if icomp>40
            disp('Reducing ICA components to 30')
            icomp = 40;
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
        close all
    end
    
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data    = ft_rejectcomponent(cfg, comps,data);
    
    
    %%
    
    cfg = [];
    cfg.continuous = 'yes';
    cfg.trl = [data.sampleinfo, 0];
    cfg.artfctdef.muscle.channel = ...
        data.label(strncmp(data.label,'MLT',3) | strncmp(data.label,'MRT',3)); % should I just look in temporal channels?
    cfg.artfctdef.muscle.cutoff  = 6;
    cfg.artfctdef.muscle.trlpadding = 0;
    cfg.artfctdef.muscle.fltpadding = 0;
    cfg.artfctdef.muscle.artpadding = 0.1;
    
    cfg.artfctdef.muscle.bpfilter    = 'yes';
    cfg.artfctdef.muscle.bpfreq      = [100 150]; % [110 140]
    cfg.artfctdef.muscle.bpfiltord   = 4; % previously 8
    cfg.artfctdef.muscle.bpfilttype  = 'but';
    cfg.artfctdef.muscle.hilbert     = 'yes';
    cfg.artfctdef.muscle.boxcar      = 0.2;
    
    [cgfz, artifact] = ft_artifact_muscle(cfg,data);
    % chans = strncmp(data.label,'MLT',3);
    % dataz = zscore(data.trial{1},[],2);
    % figure; plot(dataz(chans,:)'+ (0:nnz(chans)-1)*10 , 'k')
    % hold on
    % for ii = 1:size(artifact,1)
    %     plot(artifact(ii,1):artifact(ii,2),...
    %         dataz(chans,artifact(ii,1):artifact(ii,2))'+ (0:nnz(chans)-1)*10,'r')
    % end
    
    %%
    filt_order = []; % default
    data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[4 8],filt_order,'but');
    
    data.trial{1} = data_filt;
    clear data_filt
    musclecorropt = 0;
    %% Beamfomer
    
    windl = 2.5*60*f;
    
    % 10 minutes max, divide into 2 minutes sections? Obtain weights from whole
    % recording?
    icacomps = length(data.cfg.component);
    
    if musclecorropt == 1
        
        timew_ar = cell(1,4);
        for t = 1:4
            timew = (t-1)*windl + (1:windl);
            if all(timew<=data.sampleinfo(2))
                timedelc = [];
                for a = 1:size(artifact,1)
                    timedel = artifact(a,1):artifact(a,2);
                    [~,inda,~ ] = intersect(timew, timedel);
                    timedelc = cat(1,timedelc,inda);
                end
                timew(timedelc) = [];
                timew_ar{t} = timew;
            end
        end
        
        timew = 1:data.sampleinfo(2);
        timedelc = [];
        for a = 1:size(artifact,1)
            timedel = artifact(a,1):artifact(a,2);
            [~,inda,~ ] = intersect(timew, timedel);
            timedelc = cat(1,timedelc,inda);
        end
        timew(timedelc) = [];
        
        C = cov(data.trial{1}(:,timew)');
        E = svd(C);
        nchans = length(data.label);
        noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
        Cr = C + 4*noiseC;
        
        windl = zeros(1,4);
        for t = 1:4
            windl(t) = size(timew_ar{t},2);
        end
        
        timewgo = find(windl>0); % empty window
        windl = min(windl(timewgo));
        
        Crw = cell(1,length(timewgo));
        for t = timewgo
            datat = data.trial{1}(:,timew_ar{t});
            datat = datat(:,1:windl);
            C = cov(datat');
            E = svd(C);
            nchans = length(data.label);
            noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
            
            Crw{t} = C + 4*noiseC; % old normalization
            % Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value
        end
        
    else
        timew_ar = cell(1,4);
        for t = 1:4
            timew = (t-1)*windl + (1:windl);
            if all(timew<=data.sampleinfo(2))
                timew_ar{t} = timew; 
            end
        end
        
        C = cov(data.trial{1}');
        E = svd(C);
        nchans = length(data.label);
        noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
        Cr = C + 4*noiseC;
        
        windl = zeros(1,4);
        for t = 1:4
            windl(t) = size(timew_ar{t},2);
        end
        
        timewgo = find(windl>0); % empty window
        windl = min(windl(timewgo));
        
        Crw = cell(1,length(timewgo));
        for t = timewgo
            datat = data.trial{1}(:,timew_ar{t});
            datat = datat(:,1:windl);
            C = cov(datat');
            E = svd(C);
            nchans = length(data.label);
            noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
            
            Crw{t} = C + 4*noiseC; % old normalization
            % Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value
        end
    end
    
    W = cell(nnz(grid.inside),1);
    
    L = grid.leadfield(grid.inside);
    for ii = 1:nnz(grid.inside)
        
        lf = L{ii}; % Unit 1Am
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        W{ii} = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        
        if mod(ii,500)
            clc; fprintf('SAM weights done %.1f perc.\n',ii/nnz(grid.inside)*100)
        end
    end
    
    for t = timewgo
        Q = zeros(nnz(grid.inside),1);
        for ii = 1:nnz(grid.inside)
            w = W{ii};
            Q(ii) = w'*Crw{t}*w;
        end
        
        switch t
            case 1
                Q1 = Q;
            case 2
                Q2 = Q;
            case 3
                Q3 = Q;
            case 4
                Q4 = Q;
        end
        clear Q
    end
    
    
    if any(timewgo == 2)
        Tstat2 = (Q2-Q1)./(Q2+Q1);
        save([processing_folder,'/SAM_theta.mat'],'Tstat2');
    end
    
    if any(timewgo == 3)
        Tstat3 = (Q3-Q1)./(Q3+Q1);
        save([processing_folder,'/SAM_theta.mat'],'Tstat3','-append');
    end
    
    if any(timewgo == 4)
        Tstat4 = (Q4-Q1)./(Q4+Q1);
        save([processing_folder,'/SAM_theta.mat'],'Tstat4','-append');
    end
    
    
    
end



%%
Tstat12 = cell(1,length(param_list));
Tstat13 = cell(1,length(param_list));
Tstat14 = cell(1,length(param_list));

for zz = 1:length(param_list)
    
    clear Tstat2 Tstat3 Tstat4 grid
    data_name = param_list{zz};
    sub = data_name(1:5);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    processing_folder = [data_path,data_name,'/beamforming'];
    
    load([processing_folder,'/SAM_gamma2.mat']); % gamma 2 has better muscle correction
    load([processing_folder,'/leadfields_5mm_MNI.mat']);
    
    if exist('Tstat2','var')
        Tstat12{zz} = zeros(size(grid.inside));
        Tstat12{zz}(grid.inside) = Tstat2;
    end
    
    if exist('Tstat3','var')
        Tstat13{zz} = zeros(size(grid.inside));
        Tstat13{zz}(grid.inside) = Tstat3;
    end
    
    if exist('Tstat4','var')
        Tstat14{zz} = zeros(size(grid.inside));
        Tstat14{zz}(grid.inside) = Tstat4;
    end
    
end

Tstat2 = cell2mat(Tstat12);
Tstat3 = cell2mat(Tstat13);
Tstat4 = cell2mat(Tstat14);

if ~exist('mri_mni','var')
    mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
    % mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii');
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    gridres = 5;
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    
    atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
    
    atlas = ft_convert_units(atlas,sourcemodel.unit);
    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = 'tissue';
    sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);
    
%     aalfront = strncmp(sourcemodelAAL.tissuelabel,'Front',5) | ...
%         strncmp(sourcemodelAAL.tissuelabel,'Cingulum_Ant',12) | ...
%         strncmp(sourcemodelAAL.tissuelabel,'Rectus',6);
%   
%     Tfront  = false(size(grid.inside));
%     for ii = find(aalfront)
%         ind = find(sourcemodelAAL.tissue == ii);
%         Tfront(ind) = true;
%         %     locsAAL{ii} = sourcemodel.pos(ind,:);
%     end
    
   
    
end



aalfront([29:32,37:38,41:42,73:78]) = true;

figure; set(gcf,'position',[1175         -15         434         905])
Tfront  = false(size(grid.inside));
jj = 0;
for ii = find(aalfront)
    jj = jj+1;
    subplot(7,2,jj)
    ind = find(sourcemodelAAL.tissue == ii);
    Tfront(ind) = true;
    x2 = mean(Tstat2(ind,:));
    x3 = mean(Tstat3(ind,:));
    x4 = mean(Tstat4(ind,:));
    hold all
    errorbar(2:4,[mean(x2),mean(x3),mean(x4)],[std(x2),std(x3),std(x4)],'-o')
    xlim([1 5]); ylim([-.35 0.1]); grid on
    title(sourcemodelAAL.tissuelabel{ii})
end
set(gcf,'name','\gamma band')

% Check other ROIs
aalfront([29:32,37:38,41:42,73:78]) = true;

aalcomp = strncmp(sourcemodelAAL.tissuelabel,'Cerebellum',10) | ...
    strncmp(sourcemodelAAL.tissuelabel,'Vermis',6);
aalcomp = strncmp(sourcemodelAAL.tissuelabel,'Temp',4) ;
ind = [];
for ii = find(aalcomp)
    ind = cat(1,find(sourcemodelAAL.tissue == ii));
end
figure;
x2 = mean(Tstat2(ind,:));
x3 = mean(Tstat3(ind,:));
x4 = mean(Tstat4(ind,:));
errorbar(2:4,[mean(x2),mean(x3),mean(x4)],[std(x2),std(x3),std(x4)],'-o')

%%
close all
for jj = 1:6
for ii = 4 %2:4
    
%     Tstat = mean(eval(sprintf('Tstat%.0f',ii)),2);
    
        Tstat = eval(sprintf('Tstat%.0f',ii));
        Tstat = Tstat(:,jj);
    
    %     Tstat =  zeros(size(grid.inside));
    %     Tstat(Tfront) = Tstat0(Tfront);
    
    sourceant =[];
    sourceant.pow = Tstat;
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    
    crang = [];
    cfg = [];
    cfg.method        = 'slice';
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    
    ft_sourceplot(cfg, sourceout_Int);
    
end
end