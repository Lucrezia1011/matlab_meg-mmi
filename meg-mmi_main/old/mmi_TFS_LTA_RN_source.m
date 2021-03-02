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

Ypn = [];
Yp = [];
Yn = [];

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
    
    answer_match = bv_match.answer;
    choice_match =  bv_match.choice;
    outcome_match  = bv_match.outcome;
    mood_match = bv_match.ratemood;
    blockmood_match = bv_match.blockmood;
    slider_match = bv_match.slider;
    blockslider_match = bv_match.blockslider;
    ITI_match = bv_match.ITI ;
    buttonpress = bv_match.buttonpress;
    tasktime = bv_match.time;
    
    
    
    %% Filter data
    filt_order = []; % default
%     data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[1 150],filt_order,'but');
    data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 5,[],'but');
   
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
       outcome_match.sample(outcome_match.win==1)], datave, tasktime, [-.5 1]);
           
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

%     subplot(142); imagesc(negnor); caxis([-1 1]*2e-13); title('NEG - Neutral')
%     subplot(143); imagesc(posnor); caxis([-1 1]*2e-13); title('POS - Neutral')
%     subplot(144); imagesc(posneg); caxis([-1 1]*2e-13); title('POS - NEG')
%     set(gcf,'color','w','name',sub)
        
  
    Ypn = cat(3,Ypn,posneg);    
    Yp = cat(3,Yp,posnor);
    Yn = cat(3,Yn,negnor);
%   
    Npn = cat(1,Npn,Nposneg);
    Np  = cat(1,Np,Nposnor);
    Nn  = cat(1,Nn,Nnegnor);

    
end


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



figure;
subplot(121)
imagesc(mean(Yall,3));
% caxis([-1 1]*2e-13)
title('Average response to all feedbacks')
subplot(122)
imagesc(mean(Ypn,3));
% caxis([-1 1]*1e-13)
title('Average Positive - Negative feedback')

time = linspace(-0.5,1,npoints);
figure;
n = 10;
plot(time, mean(Yall(n,:,:),3));
hold on
plot(time, mean(Ypn(n,:,:),3));

title(sourcemodelAAL.tissuelabel{n})

save('/data/MBDU/MEG_MMI3/results/mmiRN_aal_mu5max_Z.mat','Ypn','Yp','Yn','Npn','Np','Nn');

%% Try running GLM separately for each latent variable, there is a lot of correlation between them!

% Check correlation of variables
glm1 = cell(nroi,npoints);
glm2 = cell(nroi,npoints);
glm3 = cell(nroi,npoints);
% glm4 = cell(nchans,npoints);
% glm5 = cell(nchans,npoints);
% glm6 = cell(nchans,npoints);
% glm7 = cell(nroi,npoints);
% glm8 = cell(nroi,npoints);

l = size(glm1);
N = nroi*npoints;

clc
parfor n = 1:N
%     t = 0.25;
%     [~,tt] = min(abs (dataout.time{1} - t));  
    [sen,tt] = ind2sub(l,n);
    
%     X = table(Y,RPEall,Sall,'VariableNames',{'MEG','RPE','subject'});
    
    X = ltn;
    X.MEG = squeeze(Yp(sen,tt,:));
    X.happy = happy';       
    % Try one variable at the time 
%     G1 = fitglm(X,'MEG ~ a');  
%     X = ltn.b;
    G1 = fitglm(X,'MEG ~ b ','Weights',sqrt(Np));    

    X.MEG = squeeze(Yn(sen,tt,:));
    G2 = fitglm(X,'MEG ~ b ','Weights',sqrt(Nn));  
    X.MEG = squeeze(Ypn(sen,tt,:));
    G3 = fitglm(X,'MEG ~ b ','Weights',sqrt(Npn));  
% F-stat for cluster statistics??
%     a = G1.devianceTest;
%     Fsts = a.FStat; 
        
    
%     glm1{n} =  [G1.Coefficients.Estimate(end)*1e15, G1.Coefficients.pValue(end)];
%     glm2{n} =  [G2.Coefficients.Estimate(end)*1e15, G2.Coefficients.pValue(end)];
    
    glm1{n} =  [G1.Coefficients.Estimate(2), G1.Coefficients.pValue(2)];
    glm2{n} =  [G2.Coefficients.Estimate(2), G2.Coefficients.pValue(2)];    
    glm3{n} =  [G3.Coefficients.Estimate(2), G3.Coefficients.pValue(2)];
%     glm4{n} =  [G4.Coefficients.Estimate(end)*1e15, G4.Coefficients.pValue(end)];
%     glm5{n} =  [G5.Coefficients.Estimate(end)*1e15, G5.Coefficients.pValue(end)];
%     glm6{n} =  [G6.Coefficients.Estimate(end)*1e15, G6.Coefficients.pValue(end)];

%     lme{tt} = fitlme(tbl,'MEG~RPE+(1|subject)'); %Fixed effects for RPE 
    

end
fprintf('Done')


%% PLot GLM

varnames = {'Positive Feedback (Z-score)';...
    'Negative Feedback (Z-score)'; 'Positive Feedback vs Negative (Z-score)'};

time = linspace(-.1,0.5,npoints);


for x = 1:3
    glm = eval(genvarname(['glm',num2str(x)]));
    G = cell2mat(glm);
    G = reshape(G,[nroi,length(glm{1}),npoints]);

    Gc = squeeze(G(:,1,:));
    Gp = squeeze(G(:,2,:));
%     Gc = squeeze(G(2,:,1,:));
%     Gp = squeeze(G(2,:,2,:));
      
    [sen, tt ]=  find(Gp<(0.05/(nroi)));
    [senu,ii,jj] = unique(sen);
    figure; set(gcf,'color','w','name',varnames{x})
for sn = 1:length(senu)
    switch x
        case 1
            Y = mean(Yp(senu(sn),:,:),3);
        case 2
            Y = mean(Yn(senu(sn),:,:),3);
        case 3
            Y = mean(Ypn(senu(sn),:,:),3);        
    end
    
    tx = tt(jj==sn);
    subplot(7,4,sn)
    plot(time,Y,'k')
    hold on 
    
%     [ss, ind ]=  sort(squeeze(Gp(x,senu(sn),:)));
%     ind(ss<(0.05/nroi/npoints*(1:npoints)'))
%     
    
    [g,zz] = max(abs(Gc(senu(sn),tx)));
   
    text(time(tx(zz))+0.05,Y(tx(zz)),num2str(Gc(senu(sn),tx(zz))))
    
%     tx = time(tx);
%     fill([tx(1) tx(end) tx(end) tx(1)],[-1 -1 1 1]*8e-14,[0 0 1],'facealpha',0.2,'edgecolor','none')
    plot(time(tx),Y(tx),'*r')
    title(sourcemodelAAL.tissuelabel{senu(sn)})
%     title(changroup(senu(sn))
end
end    