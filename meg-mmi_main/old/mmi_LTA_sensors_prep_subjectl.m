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
Ym = [];

Ypn = [];
Yp = [];
Yn = [];

Npn = [];
Np  = [];
Nn  = [];

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
    
   
    %%
 
    [dataout1,ttdel]= define_trials([outcome_match.sample(outcome_match.win==-1),...
       outcome_match.sample(outcome_match.win==1)], data, tasktime, [-3 3]);
           
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
    
    nroi = length(data.label);
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

%% Write latent variables

if ~exist('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/latent_vars.csv','file')
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
end

figure;
subplot(131)
imagesc(mean(Yp,3));
% caxis([-1 1]/2)
title('Average Positive feedback')
subplot(132)
imagesc(mean(Yn,3));
% caxis([-1 1]/2)
title('Average Negative feedback')
subplot(133)
imagesc(mean(Ypn,3));
% caxis([-1 1]/2)
title('Average Positive - Negative feedback')


% save('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/mmiRN_aal_mu5max_Z.mat','Ypn','Yp','Yn','Npn','Np','Nn');
Yp = reshape(Yp, [nroi*npoints, size(Yp,3)] ); % nrois x npoints
Yn = reshape(Yn, [nroi*npoints, size(Yn,3)] ); % nrois x npoints
Ypn = reshape(Ypn, [nroi*npoints, size(Ypn,3)] ); % nrois x npoints

dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/positive_sensors.txt',Yp)
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/negative_sensors.txt',Yn)
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/posneg_sensors.txt',Ypn)

