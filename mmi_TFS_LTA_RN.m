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

Ypn = [];
Yp = [];
Yn = [];

Ypl_pn = [];
Ypl_n  =[];
Ypl_p  =[];


Npn = [];
Np  = [];
Nn  = [];

%% Find common channels
sn = 2;
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds'];
h = ft_read_sens(data_name,'senstype','MEG');
label2 = h.label(strncmp(h.label,'M',1));

sn = 4;
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds'];
h = ft_read_sens(data_name,'senstype','MEG');
label4 = h.label(strncmp(h.label,'M',1));

channels = intersect(label2,label4);

%%

for sn = [1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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


    %% Clean data with ICA
    
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = channels;
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
    data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 20,[],'but');
   
    data.trial{1} = data_filt;
    clear data_filt    
    
   [dataout1,ttdel]= define_trials([outcome_match.sample(outcome_match.win==-1),...
       outcome_match.sample(outcome_match.win==1)], data, tasktime, [-.1 0.5]);
           
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
  
        dataout.trial = cat(2,dataout.trial,dataout1.trial);
        dataout.time = cat(2,dataout.time,dataout1.time);      
    end
    clear dataout1
    
    
    end
 
    idx = kmeans(X',3,'Replicates',30);
    
    figure; subplot(141)
    scatter(X(1,idx==1),X(2,idx==1))
    hold all
    scatter(X(1,idx==2),X(2,idx==2))
    scatter(X(1,idx==3),X(2,idx==3))
    xlabel('win/loss'); ylabel('RPE')
    
    [~,idxs] = sort([sum(mean(X(:,idx==1),2)), sum(mean(X(:,idx==2),2)),  sum(mean(X(:,idx==3),2))]);
    
    trials_neg = idx==idxs(1);
    trials_nor = idx==idxs(2);
    trials_pos = idx==idxs(3);
    
    nchans = length(data.label);
    
    % Add covariate noise from averaging. 
    data_neg = dataout.trial(trials_neg);
    data_neg = cell2mat(data_neg);
    data_neg = reshape(data_neg,nchans,npoints,nnz(trials_neg));
    
    data_nor = dataout.trial(trials_nor);
    data_nor = cell2mat(data_nor);
    data_nor = reshape(data_nor,nchans,npoints,nnz(trials_nor));

    data_pos = dataout.trial(trials_pos);
    data_pos = cell2mat(data_pos);
    data_pos = reshape(data_pos,nchans,npoints,nnz(trials_pos));
    
        
    avgFIC = ft_timelockanalysis([],dataout);
  
    cfg                 = [];
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, avgFIC);
    cfg.planarmethod    = 'sincos';

    %%%%%%%%%%%%%%%
%     avgFIC.avg = posneg;
%     
%     cfg.planarmethod    = 'sincos';
%     avgFICplanar        = ft_megplanar(cfg, avgFIC);
%     avgFICplanarComb = ft_combineplanar([],avgFICplanar);
%     
%     Ypl_pn = cat(3,Ypl_pn,avgFICplanarComb.avg);
    %%%%%%%%%%%%%%%
%     avgFIC.avg = posnor;
%     
%     avgFICplanar        = ft_megplanar(cfg, avgFIC);
%     avgFICplanarComb = ft_combineplanar([],avgFICplanar);
%     
%     Ypl_p = cat(3,Ypl_p,avgFICplanarComb.avg);
%     %%%%%%%%%%%%%%%%
%     avgFIC.avg = negnor;
%     
%     avgFICplanar        = ft_megplanar(cfg, avgFIC);
%     avgFICplanarComb = ft_combineplanar([],avgFICplanar);
%     
%     Ypl_n = cat(3,Ypl_n,avgFICplanarComb.avg);
    
        
    if nnz(trials_neg) < nnz(trials_pos)
        [~,ind] = sort(X(2,trials_pos),'descend');
        posneg =  mean(data_pos(:,:,ind(1:nnz(trials_neg))),3) - mean(data_neg,3);
        Nposneg = nnz(trials_neg);
        
        avgFIC.avg = mean(data_pos(:,:,ind(1:nnz(trials_neg))),3);  
        avgFICplanar        = ft_megplanar(cfg, avgFIC);
        avgFICplanarpos = ft_combineplanar([],avgFICplanar);
    
        avgFIC.avg = mean(data_neg,3);  
        avgFICplanar        = ft_megplanar(cfg, avgFIC);
        avgFICplanarneg = ft_combineplanar([],avgFICplanar);
        
        Ypl_pn = cat(3,Ypl_pn,avgFICplanarpos.avg - avgFICplanarneg.avg);
        
    elseif nnz(trials_neg) > nnz(trials_pos)
        [~,ind] = sort(X(2,trials_neg),'ascend');
        posneg = mean(data_pos,3) - mean(data_neg(:,:,ind(1:nnz(trials_pos))),3) ;
        Nposneg = nnz(trials_pos);
        
        avgFIC.avg = mean(data_pos,3);  
        avgFICplanar        = ft_megplanar(cfg, avgFIC);
        avgFICplanarpos = ft_combineplanar([],avgFICplanar);
    
        avgFIC.avg = mean(data_neg(:,:,ind(1:nnz(trials_pos))),3);  
        avgFICplanar        = ft_megplanar(cfg, avgFIC);
        avgFICplanarneg = ft_combineplanar([],avgFICplanar);
        
        Ypl_pn = cat(3,Ypl_pn,avgFICplanarpos.avg - avgFICplanarneg.avg);
        
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
     
%     subplot(142); imagesc(negnor); caxis([-1 1]*2e-13); title('NEG - Neutral')
%     subplot(143); imagesc(posnor); caxis([-1 1]*2e-13); title('POS - Neutral')
%     subplot(144); imagesc(posneg); caxis([-1 1]*2e-13); title('POS - NEG')
%     set(gcf,'color','w','name',sub)
        
      Ypn = cat(3,Ypn,posneg);
%     Yp = cat(3,Yp,posnor);
%     Yn = cat(3,Yn,negnor);
%   
    Npn = cat(1,Npn,Nposneg);
%     Np  = cat(1,Np,Nposnor);
%     Nn  = cat(1,Nn,Nnegnor);

    
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

return
w  = zeros(1,1,11);
w(1,1,:) = sqrt(Npn);

figure;
subplot(221)
avgFIC.avg = mean(Ypn,3);
cfg = [];
cfg.xlim = [0.18 0.23];
cfg.zlim = [-1 1]*4e-14;
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotER(cfg,avgFIC)

subplot(222)
avgFIC.avg = mean(Ypn.*w,3)/mean(w,3);
ft_topoplotER(cfg,avgFIC)

subplot(223)
cfg.zlim = [-1 1]*2e-14;
avgFIC.avg = mean(Ypl_pn,3);
ft_topoplotER(cfg,avgFIC)

subplot(224)
avgFIC.avg = mean(Ypl_pn.*w,3)/mean(w,3);
ft_topoplotER(cfg,avgFIC)


%% Try running GLM separately for each latent variable, there is a lot of correlation between them!

% Check correlation of variables
glm1 = cell(nchans,npoints);
glm2 = cell(nchans,npoints);
% glm3 = cell(nchans,npoints);
% glm4 = cell(nchans,npoints);
% glm5 = cell(nchans,npoints);
% glm6 = cell(nchans,npoints);
% glm7 = cell(nroi,npoints);
% glm8 = cell(nroi,npoints);


l = size(glm1);
N = nchans*npoints;

clc
parfor n = 1:N
%     t = 0.25;
%     [~,tt] = min(abs (dataout.time{1} - t));  
    [sen,tt] = ind2sub(l,n);
    
%     X = table(Y,RPEall,Sall,'VariableNames',{'MEG','RPE','subject'});
    
    X = ltn;
    X.MEG = squeeze(Ypn(sen,tt,:));
    X.happy = happy';       
    % Try one variable at the time 
%     G1 = fitglm(X,'MEG ~ a');  
%     X = ltn.b;
    G1 = fitglm(X,'MEG ~ w_RPE ','Weights',sqrt(Npn));
    
    X.MEG = squeeze(Ypl_pn(sen,tt,:));
    G2 = fitglm(X,'MEG ~ w_RPE ','Weights',sqrt(Npn));
    

% F-stat for cluster statistics??
%     a = G1.devianceTest;
%     Fsts = a.FStat; 
        
    
%     glm1{n} =  [G1.Coefficients.Estimate(end)*1e15, G1.Coefficients.pValue(end)];
%     glm2{n} =  [G2.Coefficients.Estimate(end)*1e15, G2.Coefficients.pValue(end)];
    
    glm1{n} =  [G1.Coefficients.Estimate(2)*1e15, G1.Coefficients.pValue(2)];
    glm2{n} =  [G2.Coefficients.Estimate(2)*1e15, G2.Coefficients.pValue(2)];
    
%     glm3{n} =  [G3.Coefficients.Estimate(end)*1e15, G3.Coefficients.pValue(end)];
%     glm4{n} =  [G4.Coefficients.Estimate(end)*1e15, G4.Coefficients.pValue(end)];
%     glm5{n} =  [G5.Coefficients.Estimate(end)*1e15, G5.Coefficients.pValue(end)];
%     glm6{n} =  [G6.Coefficients.Estimate(end)*1e15, G6.Coefficients.pValue(end)];

%     lme{tt} = fitlme(tbl,'MEG~RPE+(1|subject)'); %Fixed effects for RPE 
    

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

varnames = {'Positive Feedback vs Negative';'Positive Feedback vs Negative, Planar';...
    'Positive Feedback vs Neutral'; 'Negative Feedback vs Neutral'; ...
    'Positive Feedback vs Neutral, Planar'; ...
   'Negative Feedback vs Neutral, Planar'};

time = linspace(-.1,0.5,npoints);


for x = 1:2
    glm = eval(genvarname(['glm',num2str(x)]));
    G = cell2mat(glm);
    G = reshape(G,[nchans,length(glm{1}),npoints]);

    Gc = squeeze(G(:,1,:));
    Gp = squeeze(G(:,2,:));
%     Gc = squeeze(G(2,:,1,:));
%     Gp = squeeze(G(2,:,2,:));
      
    [sen, tt ]=  find(Gp<(0.05/(nchans)));
    [senu,ii,jj] = unique(sen);
    figure; set(gcf,'color','w','name',varnames{x})
for sn = 1:length(senu)
    switch x
        case 1
            Y = mean(Ypn(senu(sn),:,:),3);
        case 2
            Y = mean(Ypl_pn(senu(sn),:,:),3);
        case 3
            Y = mean(Ypn(senu(sn),:,:),3);
        case 4
            Y = mean(Ypl_p(senu(sn),:,:),3);
        case 5
            Y = mean(Ypl_n(senu(sn),:,:),3);
        case 6
            Y = mean(Ypl_pn(senu(sn),:,:),3);
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
    title(channels(senu(sn)))
%     title(changroup(senu(sn))
end
end    