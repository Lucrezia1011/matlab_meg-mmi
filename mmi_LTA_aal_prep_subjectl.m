clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

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

Yall = [];
Ym = [];

Ypn = [];
Yp = [];
Yn = [];

Ypns = [];
Yps = [];
Yns = [];

Npn = [];
Np  = [];
Nn  = [];

%%

for sn =[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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
    
    X = [];
    dataout = [];
    for runs = 1:length(data_names)
    
    %%
    data_name = data_names{runs};
    
    sub = data_name(1:5);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
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
    locsAAL = zeros(R,3);
    for ii = 1:R
        ind = find(sourcemodelAAL.tissue == ii);
        locs(ii,:) = mean(sourcemodel.pos(ind,:));
        
        locsAAL(ii,:) = mean(template_grid.pos(ind,:));
    end
    
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
        
    end
    
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data          = ft_rejectcomponent(cfg, comps,data);
    
      %% Read events
    
%     bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
%     for ii = 1:length(bv_names)
%         if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
%             bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
%         end
%     end
           
    bv_match = match_triggers_fc(data_name);
    
%     answer_match = bv_match.answer;
%     choice_match =  bv_match.choice;
    outcome_match  = bv_match.outcome;
%     mood_match = bv_match.ratemood;
%     blockmood_match = bv_match.blockmood;
%     slider_match = bv_match.slider;
%     blockslider_match = bv_match.blockslider;
%     ITI_match = bv_match.ITI ;
%     buttonpress = bv_match.buttonpress;
    tasktime = bv_match.time;
    
    
    
    %% Filter data
    filt_order = []; % default
%     data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[1 150],filt_order,'but');
    data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35,[],'but');
   
    data.trial{1} = data_filt;
    clear data_filt    
    
    %% Beamfomer
    icacomps = length(data.cfg.component);
    
    C = cov(data.trial{1}');
    E = svd(C);
    nchans = length(data.label);
    noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
    
    Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value
%     Cr = C + 4*noiseC; % x4 min singular value
 
    L = grid.leadfield;
    clear VE
    
    VE(1:size(L,2)) = {0};
    W(1:size(L,2)) = {0};
    
    
    for ii = 1:length(L)
        lf = L{ii}; % Unit 1Am
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
       
        lfo = lf*v(:,jj); % Lead field with selected orientation
       
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        W{ii} = w;
        VE{ii}  = w'*data.trial{1};
        
        clc
        fprintf('SAM running %d/%d .\n', ii, R)
        
    end
    clc
    fprintf('Beamformer finished\n' )
    
    VE = cell2mat(VE');
    VE = zscore(VE,0,2);
    %%
    datave = data;
    datave.trial{1} = VE;
    datave.label = sourcemodelAAL.tissuelabel';
    
    [dataout1,ttdel]= define_trials([outcome_match.sample(outcome_match.win==-1),...
       outcome_match.sample(outcome_match.win==1)], datave, tasktime, [-3 3]);
           
    cfg = [];
    cfg.resamplefs = 200; % Downsample to 200Hz for ease of memory
    dataout1 = ft_resampledata(cfg, dataout1);
    
    npoints = length(dataout1.time{1});
    
    % Separate outcomes into loss, win and neutral
    x = [outcome_match.loseamount(outcome_match.win==-1),outcome_match.winamount(outcome_match.win==1);...
        outcome_match.RPE(outcome_match.win==-1) ,outcome_match.RPE(outcome_match.win==1)];
    x(:,ttdel) = [];
    
    X = cat(2,X,x);
    if isempty(dataout)
        dataout = dataout1;       
    else 
        avgFIC1 = ft_timelockanalysis([],dataout1);
        avgFIC  = ft_timelockanalysis([],dataout);      
        vec = corr(avgFIC.avg',avgFIC1.avg');
        for ii = 1:length(dataout1.trial)
            dataout1.trial{ii} = dataout1.trial{ii}.*sign(diag(vec));            
        end
        
        dataout.trial = cat(2,dataout.trial,dataout1.trial);
        dataout.time = cat(2,dataout.time,dataout1.time);      
    end
    clear dataout1
    
    
    end
 
    idx = kmeans(X',3,'Replicates',30);
    
%     figure; subplot(141)
%     scatter(X(1,idx==1),X(2,idx==1))
%     hold all
%     scatter(X(1,idx==2),X(2,idx==2))
%     scatter(X(1,idx==3),X(2,idx==3))
%     xlabel('win/loss'); ylabel('RPE')
    
    [~,idxs] = sort([sum(mean(X(:,idx==1),2)), sum(mean(X(:,idx==2),2)),  sum(mean(X(:,idx==3),2))]);
    
    trials_neg = idx==idxs(1);
    trials_nor = idx==idxs(2);
    trials_pos = idx==idxs(3);
    
    nroi = length(datave.label);
    % Add covariate noise from averaging. 
    data_neg = dataout.trial(trials_neg);
    data_neg = cell2mat(data_neg);
    data_neg = reshape(data_neg,nroi,npoints,nnz(trials_neg));
    
    data_nor = dataout.trial(trials_nor);
    data_nor = cell2mat(data_nor);
    data_nor = reshape(data_nor,nroi,npoints,nnz(trials_nor));

    data_pos = dataout.trial(trials_pos);
    data_pos = cell2mat(data_pos);
    data_pos = reshape(data_pos,nroi,npoints,nnz(trials_pos));
    

    Y = cat(3,data_neg,data_nor,data_pos);
    Ys = Y;    
   % for evoked responses!! To deal with sign uncertainty    
    if ~isempty(Yall)
        vec = corr(mean(Y,3)');
        Y = Y.*sign(vec(:,69)); % aligns to the left precentral lobule
        vec = corr(mean(Y,3)',mean(Yall,3)');
        Y = Y.*sign(vec(69,69));
    else
        vec = corr(mean(Y,3)');
        Y = Y.*sign(vec(:,69));
    end
    Yall = cat(3,Yall,Y);
       
    Ym = cat(3,Ym,mean(Ys,3));
   
    
    data_neg = Y(:,:,1:nnz(trials_neg));
    data_pos = Y(:,:,(end-nnz(trials_pos)+1):end);
    data_nor = Y(:,:,(nnz(trials_neg)+1):(end-nnz(trials_pos)));
    
    if nnz(trials_neg) < nnz(trials_pos)
        [~,ind] = sort(X(2,trials_pos),'descend');
        posneg =  mean(data_pos(:,:,ind(1:nnz(trials_neg))),3) - mean(data_neg,3);
        Nposneg = nnz(trials_neg);
      
    elseif nnz(trials_neg) > nnz(trials_pos)
        [~,ind] = sort(X(2,trials_neg),'ascend');
        posneg = mean(data_pos,3) - mean(data_neg(:,:,ind(1:nnz(trials_pos))),3) ;
        Nposneg = nnz(trials_pos);
       
    end
    
    data_negs = Ys(:,:,1:nnz(trials_neg));
    data_poss = Ys(:,:,(end-nnz(trials_pos)+1):end);
    
    if nnz(trials_neg) < nnz(trials_pos)
        [~,ind] = sort(X(2,trials_pos),'descend');
        posnegs =  mean(data_poss(:,:,ind(1:nnz(trials_neg))),3) - mean(data_negs,3);
        Nposneg = nnz(trials_neg);
      
    elseif nnz(trials_neg) > nnz(trials_pos)
        [~,ind] = sort(X(2,trials_neg),'ascend');
        posnegs = mean(data_poss,3) - mean(data_negs(:,:,ind(1:nnz(trials_pos))),3) ;
        Nposneg = nnz(trials_pos);
       
    end
    
    
%     if nnz(trials_neg) < nnz(trials_nor)
%         [~,ind] = sort(abs(X(2,trials_nor)),'ascend');
%         negnor = mean(data_neg,3) - mean(data_nor(:,:,ind(1:nnz(trials_neg))),3);
%         Nnegnor = nnz(trials_neg);
%     elseif nnz(trials_neg) > nnz(trials_nor)
%         [~,ind] = sort(X(2,trials_neg),'ascend');
%         negnor = mean(data_neg(:,:,ind(1:nnz(trials_nor))),3) - mean(data_nor,3);
%         Nnegnor = nnz(trials_nor);
%     end
%     
%     if nnz(trials_pos) < nnz(trials_nor)
%         [~,ind] = sort(abs(X(2,trials_nor)),'ascend');
%         posnor = mean(data_pos,3) - mean(data_nor(:,:,ind(1:nnz(trials_pos))),3);
%         Nposnor = nnz(trials_pos);
%     elseif nnz(trials_pos) > nnz(trials_nor)
%         [~,ind] = sort(X(2,trials_pos),'descend');
%         posnor = mean(data_pos(:,:,ind(1:nnz(trials_nor))),3) - mean(data_nor,3);
%         Nposnor = nnz(trials_nor);
%     end
%      
     posnor = mean(data_pos,3);
     Nposnor = nnz(trials_pos);

     negnor = mean(data_neg,3);
     Nnegnor = nnz(trials_neg);

     posnors = mean(data_poss,3);
     negnors = mean(data_negs,3);
     
     
%     subplot(142); imagesc(negnor); caxis([-1 1]*2e-13); title('NEG - Neutral')
%     subplot(143); imagesc(posnor); caxis([-1 1]*2e-13); title('POS - Neutral')
%     subplot(144); imagesc(posneg); caxis([-1 1]*2e-13); title('POS - NEG')
%     set(gcf,'color','w','name',sub)
        
  
    Ypn = cat(3,Ypn,posneg);    
    Yp = cat(3,Yp,posnor);
    Yn = cat(3,Yn,negnor);
%   

    Ypns = cat(3,Ypns,posnegs);    
    Yps = cat(3,Yps,posnors);
    Yns = cat(3,Yns,negnors);



    Npn = cat(1,Npn,Nposneg);
    Np  = cat(1,Np,Nposnor);
    Nn  = cat(1,Nn,Nnegnor);

    
end


%% Check difference between sign flip methods!

S = zeros(size(Ym));
npoints = size(Ym,2);
for ii = 1:size(Ym,1)
    Yu = squeeze(Ym(ii,:,:))';

    %    % first pass
    Yav = [];
    C = corr(Yu');
    C = triu(C,1);
    C( abs(C-1) < 0.01) = 0;
    [~,ind] = sort(abs(C(:)),'descend');
    m = sign(C(ind(1)));
    [i,j] = ind2sub(size(C),ind(1));
    Yav = cat(1,Yav, Yu(i,:), m*Yu(j,:));
    s = zeros(size(C,1),1);
    s(i) = 1;
    s(j) = m;
    
    while nnz(s) < size(Yu,1)
        C = corr(mean(Yav,1)',Yu');
        [~,ind] = sort(abs(C(:)),'descend');
        inds = find(s);
        z = 1;
        while any(ind(z) == inds)
            z = z+1;
        end
        
        m = sign(C(ind(z)));
        s(ind(z)) = m;
        Yav =  cat(1,Yav, m*Yu(ind(z),:));
    end
    
%     time = linspace(-0.5,1,size(Yav,2));
%     figure; plot(time,mean(Yav,1))

    S(ii,:,:) = repmat(s',[npoints,1]);
  
end

Ym = Ym.*S;
C = corr(mean(Ym,3)');
s = sign(C(:,69));
S = S.*repmat(s,[1,npoints,size(S,3)]);

Ypns = Ypns.*S;
Yps = Yps.*S;
Yns = Yns.*S;


%%
bv_name = '/data/MBDU/MEG_MMI3/data/behavioral/Best_Fit_Parameters_MEG.csv';

opts = detectImportOptions(bv_name);
ltv = readtable(bv_name,opts); % bahavioral data

bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');

ltn  = table;
jj = 0;
happy = zeros(1,11);
for sn = [1:7,11,14:16]
    jj = jj+1;
    ltn(jj,:) = ltv(find(ltv.subject_ID == str2num(subn(sn,:))),:);
    
    for ii = 1:length(bv_names)
        if strncmp(bv_names(ii).name,subn(sn,:),5)
            bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
        end
    end
    opts = detectImportOptions(bv_name);
    bv = readtable(bv_name,opts);
    happy(jj) = bv.lifeHappySlider_response(4); % bahavioral data
end
ltn.happy = happy';
ltn.Npos = Np;
ltn.Nneg = Nn;
ltn.Nposneg = Npn;


writetable(ltn,'/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/latent_vars.csv');


figure;
subplot(231)
imagesc(mean(Yp,3));
caxis([-1 1]/2)
title('Average Positive feedback')
subplot(232)
imagesc(mean(Yn,3));
caxis([-1 1]/2)
title('Average Negative feedback')
subplot(233)
imagesc(mean(Ypn,3));
caxis([-1 1]/2)
title('Average Positive - Negative feedback')

subplot(234)
imagesc(mean(Yps,3));
caxis([-1 1]/2)
title('Average Positive feedback')
subplot(235)
imagesc(mean(Yns,3));
caxis([-1 1]/2)
title('Average Negative feedback')
subplot(236)
imagesc(mean(Ypns,3));
caxis([-1 1]/2)
title('Average Positive - Negative feedback')


time = linspace(-3,3,npoints);
n = 5;
figure(2); clf 
subplot(211); plot(time,squeeze(Yp(n,:,:)))
hold on
plot(time,mean(squeeze(Yp(n,:,:)),2),'k','linewidth',2)
ylim([-1.5 1.5])

subplot(212); plot(time,squeeze(Yps(n,:,:)))
hold on
plot(time,mean(squeeze(Yps(n,:,:)),2),'k','linewidth',2)
ylim([-1.5 1.5])

% save('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/mmiRN_aal_mu5max_Z.mat','Ypn','Yp','Yn','Npn','Np','Nn');
Yps = reshape(Yps, [nroi*npoints, size(Ym,3)] ); % nrois x npoints
Yns = reshape(Yns, [nroi*npoints, size(Ym,3)] ); % nrois x npoints
Ypns = reshape(Ypns, [nroi*npoints, size(Ym,3)] ); % nrois x npoints

dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/positive_aal_mu5max_Z.txt',Yps)
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/negative_aal_mu5max_Z.txt',Yns)
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/posneg_aal_mu5max_Z.txt',Ypns)

save('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/mmi_aal_mu5max_Z.mat','Ypns','Yps','Yns','Npn','Np','Nn');

