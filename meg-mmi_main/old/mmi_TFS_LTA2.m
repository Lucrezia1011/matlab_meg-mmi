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

npoints = 300;
Yall = [];
Ybeta = [];
Ytheta  =[];
Sall = [];
sensall = [];
RPEall = [];
ltvall = [];

% ltvall = table([],[],[],[],[],[],[],[],[],[],[],'VariableNames',...
%         {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    


for sn = 7 %sn = [1:6,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
    
    
    clearvars -except sn subn locs Yall Ybeta Ytheta Sall sensall RPEall ltvall channels
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
    
    bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
    for ii = 1:length(bv_names)
        if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
            bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
        end
    end
    
    opts = detectImportOptions(bv_name);
    ltv = readtable(bv_name,opts); % bahavioral data
        
    
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
%     fprintf('Beamformer finished\n' )
    
    VE = cell2mat(VE');
    
    
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
    buttonpress = bv_match.buttonpress;
    tasktime = bv_match.time;
    
    
    %% Read behavioral file
    
    indbm = blockmood_match.sample~=0;
    indm = mood_match.sample~=0;
    
    [x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
    v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
    Fmood = griddedInterpolant(x,v(ind),'pchip');
    
    xi = outcome_match.sample(outcome_match.win~=0);
    moodi = Fmood(xi); % interpolated mood timecourse
    
    
    %%
    datave = data;
    datave.trial{1} = VE;
    datave.label = sourcemodelAAL.tissuelabel';
    dataf = datave;
    
    % Evoked responses 
    data_filt = ft_preproc_lowpassfilter(datave.trial{1}, datave.fsample, 35,[],'but');
    dataf.trial{1}= data_filt;
    
    [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
           
    ltvind = outcome_match.bv_index(outcome_match.win~=0) - 12; % indeces start at 13
    ltvind(ttdel) = [];
    
    cfg = [];
    cfg.resamplefs = 200; % Downsample to 300Hz for ease of memory
    dataout = ft_resampledata(cfg, dataout);
        
    
    nchans = length(dataout.label);
 
    % Try this with all subjects!
    ntrials = length(dataout.trial);
    % S = sn*ones(nchans*ntrials,1);
    S = sn*ones(ntrials,1);
    % sen = repmat(dataout.label,[ntrials,1]);
    Sall = cat(1,Sall,S);
     
    ltvcut = [ltv.EC, ltv.EG, ltv.Ediff, ltv.LTA, ltv.new_p, ltv.RPE, ltv.LTA_sum, ltv.RPE_sum];
    ltvcut = ltvcut(ltvind,:);
    trials = ltv.Var1;
    trials = trials(ltvind);
    
    ltvcut(:,4:11) = ltvcut;
    ltvcut(:,1) = S;
    ltvcut(:,2) = trials;
    ltvcut(:,3) = moodi;

    
    Y =  cell2mat(dataout.trial);
    Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);
    
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
    %%%%%%%%%%%%%%
    
    Yall = cat(3,Yall,Y);
    
     % Induced oscillations
    data_filt = ft_preproc_bandpassfilter(datave.trial{1}, datave.fsample, [4 8],[],'but');
    data_filt = abs(hilbert(data_filt'));
    dataf.trial{1}= zscore(data_filt)';
    [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
    cfg = [];
    cfg.resamplefs = 100; % Downsample to 300Hz for ease of memory
    dataout = ft_resampledata(cfg, dataout);    
    Y =  cell2mat(dataout.trial);
    Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);
    Ytheta = cat(3,Ytheta,Y);
    
    data_filt = ft_preproc_bandpassfilter(datave.trial{1}, datave.fsample, [13 30],[],'but');
    data_filt = abs(hilbert(data_filt'));
    dataf.trial{1}= zscore(data_filt)';
    [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
    cfg = [];
    cfg.resamplefs = 100; % Downsample to 300Hz for ease of memory
    dataout = ft_resampledata(cfg, dataout);    
    Y =  cell2mat(dataout.trial);
    Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);
    Ybeta = cat(3,Ybeta,Y);
 
    
    ltvall = cat(1,ltvall,ltvcut);
    end  
end

return

%% Try running GLM separately for each latent variable, there is a lot of correlation between them!

nroi = length(sourcemodelAAL.tissuelabel);

npoints = size(Ybeta,2);
glmp = cell(nroi,npoints);
glmc = cell(nroi,npoints);
glm1 = cell(nroi,npoints);
glm2 = cell(nroi,npoints);
glm3 = cell(nroi,npoints);
glm4 = cell(nroi,npoints);
glm5 = cell(nroi,npoints);
glm6 = cell(nroi,npoints);
glm7 = cell(nroi,npoints);
glm8 = cell(nroi,npoints);



l = size(glmp);
N = nroi*npoints;

clc
parfor n = 1:N
%     t = 0.25;
%     [~,tt] = min(abs (dataout.time{1} - t));  
    [sen,tt] = ind2sub(l,n);

    Y = squeeze(Ytheta(sen,tt,:));
    
%     X = table(Y,RPEall,Sall,'VariableNames',{'MEG','RPE','subject'});
    
    X = table(Y,Sall,ltvall(:,1),ltvall(:,2),ltvall(:,3),ltvall(:,4),...
        ltvall(:,5),ltvall(:,6),ltvall(:,7),ltvall(:,8),'VariableNames',...
        {'MEG','subject','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    
    % Try one variable at the time 
    G1 = fitglm(X,'MEG ~ EC + subject','CategoricalVars','subject');
    G2 = fitglm(X,'MEG ~ EG + subject','CategoricalVars','subject');
    G3 = fitglm(X,'MEG ~ Ediff + subject','CategoricalVars','subject');
    G4 = fitglm(X,'MEG ~ LTA + subject','CategoricalVars','subject');
    G5 = fitglm(X,'MEG ~ new_p + subject','CategoricalVars','subject');
    G6 = fitglm(X,'MEG ~ RPE + subject','CategoricalVars','subject');
    G7 = fitglm(X,'MEG ~ LTA_sum + subject','CategoricalVars','subject');
    G8 = fitglm(X,'MEG ~ RPE_sum + subject','CategoricalVars','subject');
    
% F-stat for cluster statistics??
%     a = G1.devianceTest;
%     Fsts = a.FStat; 
    
    G = fitlme(X,'MEG ~ LTA_sum + (LTA_sum|subject) + 1|trial');    
    
    
    glmp{n} = G.Coefficients.pValue((end-3):end);
    glmc{n} = G.Coefficients.Estimate((end-3):end);
    
    glm1{n} =  [G1.Coefficients.Estimate(end), G1.Coefficients.pValue(end)];
    glm2{n} =  [G2.Coefficients.Estimate(end), G2.Coefficients.pValue(end)];
    glm3{n} =  [G3.Coefficients.Estimate(end), G3.Coefficients.pValue(end)];
    glm4{n} =  [G4.Coefficients.Estimate(end), G4.Coefficients.pValue(end)];
    glm5{n} =  [G5.Coefficients.Estimate(end), G5.Coefficients.pValue(end)];
    glm6{n} =  [G6.Coefficients.Estimate(end), G6.Coefficients.pValue(end)];
    glm7{n} =  [G7.Coefficients.Estimate(end), G7.Coefficients.pValue(end)];
    glm8{n} =  [G8.Coefficients.Estimate(end), G8.Coefficients.pValue(end)];

%     lme{tt} = fitlme(tbl,'MEG~RPE+(1|subject)'); %Fixed effects for RPE 
    
%     fprintf('Done %d/%d\n',n,N)
end
fprintf('Done')

%% PLot GLM
Xo = table(squeeze(Ybeta(1,1,:)),Sall,ltvall(:,1),ltvall(:,2),ltvall(:,3),ltvall(:,4),...
        ltvall(:,5),ltvall(:,6),ltvall(:,7),ltvall(:,8),'VariableNames',...
        {'MEG','subject','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});


time = linspace(-0.5,1,npoints);
Gp = cell2mat(glmp);
Gp = reshape(Gp,[length(glmp{1}),size(glmp)]);

Gc = cell2mat(glmc);
Gc = reshape(Gc,[length(glmc{1}),size(glmc)]);

% varnames = Xo.Properties.VariableNames([4,7:9]);
varnames = Xo.Properties.VariableNames([4:10]);


for x = 1:7
[sen, tt ]=  find(squeeze(Gp(x,:,:))<(0.05/(npoints)));
[senu,ii,jj] = unique(sen);
figure; set(gcf,'color','w','name',varnames{x})
for sn = 1:length(senu)
%     Y = mean(Yall(senu(sn),:,:),3);
    Y = mean(Ytheta(senu(sn),:,:),3);
    tx = tt(jj==sn);
    subplot(7,4,sn)
    plot(time,Y,'k')
    hold on 
    
%     [ss, ind ]=  sort(squeeze(Gp(x,senu(sn),:)));
%     ind(ss<(0.05/nroi/npoints*(1:npoints)'))
%     
    
    [g,zz] = max(abs(Gc(x,senu(sn),tx)));
   
    text(time(tx(zz))+0.05,Y(tx(zz)),num2str(Gc(x,senu(sn),tx(zz))))
    
%     tx = time(tx);
%     fill([tx(1) tx(end) tx(end) tx(1)],[-1 -1 1 1]*8e-14,[0 0 1],'facealpha',0.2,'edgecolor','none')
    plot(time(tx),Y(tx),'*r')
    title(sourcemodelAAL.tissuelabel(senu(sn)))
%     title(changroup(senu(sn))
end
end    


%% PLot GLM
Xo = table(squeeze(Ybeta(1,1,:)),Sall,ltvall(:,1),ltvall(:,2),ltvall(:,3),ltvall(:,4),...
        ltvall(:,5),ltvall(:,6),ltvall(:,7),ltvall(:,8),'VariableNames',...
        {'MEG','subject','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});


time = linspace(-0.5,1,npoints);



varnames = Xo.Properties.VariableNames(3:end);


for x = 1:8
    glm = eval(genvarname(['glm',num2str(x)]));
    G = cell2mat(glm);
    G = reshape(G,[nroi,length(glm{1}),npoints]);

    Gc = squeeze(G(:,1,:));
    Gp = squeeze(G(:,2,:));
    
      
    [sen, tt ]=  find(Gp<(0.05/(npoints)));
    [senu,ii,jj] = unique(sen);
    figure; set(gcf,'color','w','name',varnames{x})
for sn = 1:length(senu)
%     Y = mean(Yall(senu(sn),:,:),3);
    Y = mean(Ytheta(senu(sn),:,:),3);
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
    title(sourcemodelAAL.tissuelabel(senu(sn)))
%     title(changroup(senu(sn))
end
end    