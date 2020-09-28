clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

for sn = [1,3,4,6,7,8,9,11,14,15,16] %Subjects showing enough variation in mood
    sub = subn(sn,:);
    
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    if ~exist(data_name,'dir')
        data_name = [sub,'MMI_mmi3_proc1.ds'];
    end
    
    if strcmp(sub, '24201') || strcmp(sub, '22695')
        data_name = [sub,'MMI_mmi3_proc2.ds'];
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
    
    
    answer_match = bv_match.answer;
    choice_match =  bv_match.choice;
    outcome_match  = bv_match.outcome;
    mood_match = bv_match.ratemood;
    blockmood_match = bv_match.blockmood;
    slider_match = bv_match.slider;
    blockslider_match = bv_match.blockslider;
    ITI_match = bv_match.ITI ;
    
    pRPE = outcome_match.win == 1 ;
    nRPE = outcome_match.win == -1 ;
    pRPE_sample = outcome_match.sample(pRPE);
    nRPE_sample = outcome_match.sample(nRPE);
    
    %% Filter data in theta band
    freqband = [4,8];
    
    filt_order = []; % default
    data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, freqband,filt_order,'but')';
    
    mu = 4;
    C = cov(data_filt);
    noise = min(svd(C)).*eye(size(C));
    Cr = C + mu*noise;
    
    
    %% Cue-P3 ERP: 200ms
    %
    % wind = 0.2*f:0.4*f;
    % % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    %
    % % win_logic = bv.RPE > 10 & bv.outcomeAmount > 10;
    % win_logic = bv.RPE > 0;
    % win_logic = win_logic(1:length(outcome_match.sample));
    % win_samples = outcome_match.sample(win_logic);
    %
    % % lose_logic = bv.RPE < -10 & bv.outcomeAmount < -10;
    % lose_logic = bv.RPE < 0;
    % lose_logic = lose_logic(1:length(outcome_match.sample));
    % lose_samples = outcome_match.sample(lose_logic);
    %
    % Cap = zeros([size(C),size(win_samples,2)]);
    % for tt= 1:size(win_samples,2)
    %     data_win = data_filt(win_samples(tt)+wind,:);
    %     Cap(:,:,ii) = cov(data_win);
    % end
    % Cap = mean(Cap,3);=
    %
    % Can = zeros([size(C),size(lose_samples,2)]);
    % for tt= 1:size(lose _samples,2)
    %     data_lose = data_filt(lose_samples(tt)+wind,:);
    %     Can(:,:,ii) = cov(data_lose);
    % end
    % Can = mean(Can,3);
    %

    
    %% Mood decision 3s
    
    ind = mood_match.sample>0;
    mood_match.mood(~ind)= NaN;
    
    v1 = mood_match.mood(ind);
    samples1 = mood_match.sample(ind);
    
    ind = blockmood_match.sample>0;
    blockmood_match.mood(~ind) = NaN;
    v0 = blockmood_match.mood(ind);
    samples0 = blockmood_match.sample(ind);
    
    v = [v1,v0];
    [s,ind] = sort(v);
    
    samples = [samples1, samples0];
    samples = samples(ind); 
       
    hs = s > median(s);
    ls = s < median(s);
    
    happy_samples = samples(hs);
    sad_samples = samples(ls);
    
    if nnz(hs)> nnz(ls)
        happy_samples(1) = [];
    elseif nnz(hs)< nnz(ls)
        sad_samples(end) = [];
    end
    
    wind = 0*f:3*f; % Exclude edges for sensory and motor overlap
    % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    
    
    Cah = zeros([size(C),size(happy_samples,2)]);
    for tt= 1:size(happy_samples,2)
        data_win = data_filt(happy_samples(tt)+wind,:);
        Cah(:,:,ii) = cov(data_win);
    end
    Cah = mean(Cah,3);
    
    Cal = zeros([size(C),size(sad_samples,2)]);
    for tt= 1:size(sad_samples,2)
        data_lose = data_filt(sad_samples(tt)+wind,:);
        Cal(:,:,ii) = cov(data_lose);
    end
    Cal = mean(Cal,3);
    
    
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
            Zstat_mood{ii} = ((Qah + Qal)/2 ) ./ n;

        end
        
        if mod(ii,300) == 0           
            fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
        end
        
    end
    clc
    fprintf('SAM finished\n' )
    %% Save result
    
    if ~exist([data_path,'results'],'dir')
        mkdir results
    end
    cd results
    
    mriname = 'ft_coreg_anat.nii';
    if ~exist(mriname,'file')
        ft_write_mri(mriname,mri,'dataformat','nifti');
    end
    
    
    cfg = [];
    cfg.parameter = 'pow';
    sourceTstat = struct;
    sourceTstat.dim = grid.dim;
    sourceTstat.inside = grid.inside;
    sourceTstat.pos = grid.pos;
    sourceTstat.method = 'average';
    
    %     if ~exist(zname,'file')
    
    sourceTstat.avg.pow =  cell2mat(Zstat_mood);
    sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
    sourcePostInt.anatomy = sourcePostInt.pow;
    
    zname =  sprintf('MoodZ_3s_PseudoT_%d-%dHz',freqband(1),freqband(2));
    zname_nii = [zname,'.nii'];
    ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
    
    
    sourceTstat.avg.pow =  cell2mat(Tstat_mood);
    sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
    sourcePostInt.anatomy = sourcePostInt.pow;
    
    % 3s before mood rating high - low mood
    zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz',freqband(1),freqband(2));
    zname_nii = [zname,'.nii'];
    ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
    
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
    [mrinorm] = ft_volumenormalise(cfg, mri);
    [sourcenorm] = ft_volumenormalise(cfg, sourcePostInt);
    
    zname =  sprintf('ft_coreg_anat_norm');
    zname_nii = [zname,'.nii'];
    ft_write_mri(zname_nii,mrinorm,'dataformat','nifti');
    
    sourcenorm.anatomy = sourcenorm.pow;
    zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz_norm',freqband(1),freqband(2));
    zname_nii = [zname,'.nii'];
    ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');
    
end