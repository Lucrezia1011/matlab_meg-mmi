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


% ntrials_pRPE = cell(1,16);
% ntrials_nRPE = cell(1,16);
% ntrials_outcome = cell(1,16);

timelock_trials = cell(1,16);
timelock_trialsc =cell(1,16);
timelock_choice = cell(1,16);
timelock_choicec = cell(1,16);
timelock_out = cell(1,16);
timelock_pRPE = cell(1,16);
timelock_nRPE = cell(1,16);

tfs_trials = cell(1,16);
tfs_trialsc =cell(1,16);
tfs_choice = cell(1,16);
tfs_choicec = cell(1,16);
tfs_out = cell(1,16);
tfs_pRPE = cell(1,16);
tfs_nRPE = cell(1,16);

RT_gamble =  cell(1,16);
RT_certain =  cell(1,16);

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
% Can do the same for ERF!
for sn = [1:9,11,12,14:16] % co-register 2 then run these! %[1,3,4,6,7,8,9,11,14,15,16]   %Subjects showing enough variation in mood
     
    data_button_all = [];
    data_outcome_all = [];
    data_trial_all = [];
    data_choice_all = [];   
        
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
    
    for runs = 1:length(data_names)
        
        cfg = [];
        cfg.dataset = data_names{runs};
        cfg.continuous = 'yes';
        cfg.channel = channels;
        data = ft_preprocessing(cfg);
        
        processing_folder = [data_path,data_names{runs},'/beamforming'];
        if ~exist(processing_folder,'dir')
            mkdir(processing_folder)
        end
        
        f = data.fsample;
        
        if exist([processing_folder,'/ICA_artifacts.mat'],'file')
            ICAartf = load([processing_folder,'/ICA_artifacts.mat']);
            comps = ICAartf.comps;
            cfg           = [];
            cfg.component = 1:length(comps.label);
            data        = ft_rejectcomponent(cfg, comps,data);
        end
                
        
        filt_order = []; % default
        % data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, freq_band(freq,:),filt_order,'but');
        data.trial{1} = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 ); % 35Hz for ERF

        %% Read events
        
        bv_match = match_triggers_fc(data_names{runs});
        
        answer_match = bv_match.answer;
        choice_match =  bv_match.choice;
        outcome_match  = bv_match.outcome;
        mood_match = bv_match.ratemood;
        blockmood_match = bv_match.blockmood;
        slider_match = bv_match.slider;
        blockslider_match = bv_match.blockslider;
        ITI_match = bv_match.ITI ;
        buttonpress = bv_match.buttonpress;
        %
%         pind = (outcome_match.loseamount<0 & outcome_match.winamount>0 & outcome_match.RPE>4 | outcome_match.RPE > 8);
%         nind = (outcome_match.loseamount<0 & outcome_match.RPE<-4 | outcome_match.RPE<-8);
%         
%         if nnz(pind) > nnz(nind)
%             [RPE,ind] = sort(outcome_match.RPE,'descend');
%             pind = ind(1:nnz(nind));
%         elseif nnz(pind) < nnz(nind)
%             [RPE,ind] = sort(outcome_match.RPE,'ascend');
%             nind = ind(1:nnz(pind));
%         end
            
%         [RPE,ind] = sort(outcome_match.RPE,'descend');
%         pind = ind(RPE>4);
%         [RPE,ind] = sort(outcome_match.RPE,'ascend');
%         nind = ind(RPE<-4);
%         
%         data_positive = define_trials(outcome_match.sample(pind),data,bv_match.time,[-3 4]);
%         data_negative = define_trials(outcome_match.sample(nind),data,bv_match.time,[-3 4]);
%         
%         npos = size(data_positive.trial,2);
%         nneg = size(data_negative.trial,2);
%         if npos > nneg
%             data_positive.trial = data_positive.trial(1:nneg);
%             data_positive.time = data_positive.time(1:nneg);
%             data_positive.sampleinfo = data_positive.sampleinfo(1:nneg,:);
%         elseif npos < nneg
%             data_negative.trial = data_negative.trial(1:npos);
%             data_negative.time = data_negative.time(1:npos);
%             data_negative.sampleinfo = data_negative.sampleinfo(1:npos,:);
%         end
%         data_button = define_trials(buttonpress,data,bv_match.time,[-0.2 1]);
        outcome_wind  = [-1 3];
        choice_wind  = [-2.5 2.5];    
        
        Rsamples = outcome_match.sample(outcome_match.win~=0);
        [data_outcome, ttdel] = define_trials(Rsamples,data,bv_match.time,outcome_wind);
        RPE = outcome_match.RPE(outcome_match.win~=0);
        RPE(ttdel) = [];
        outcome = outcome_match.winamount;
        outcome(outcome_match.win==-1) = outcome_match.loseamount(outcome_match.win==-1);
        outcome = outcome(outcome_match.win ~= 0);
        outcome(ttdel) = [];
        % Select trials where participant choose an option
        
        Rsamples = choice_match.sample(choice_match.sample~=0 & ~isnan(choice_match.choice));
        [data_choice,ttdel] = define_trials(Rsamples,data,bv_match.time,choice_wind);
        
        gamble_choice =  choice_match.choice(choice_match.sample~=0 & ~isnan(choice_match.choice)) == 1;
        if ~isempty(ttdel)
            gamble_choice(ttdel) = [];
        end
        
        Rsamples = answer_match.sample(answer_match.sample~=0 & ~isnan(answer_match.choice));
        [data_trial,ttdel] = define_trials(Rsamples,data,bv_match.time,outcome_wind);
        
        gamble_trial = answer_match.choice(answer_match.sample~=0 & ~isnan(answer_match.choice)) == 1;
        Rtime = answer_match.RT(answer_match.sample~=0 & ~isnan(answer_match.choice));
        if ~isempty(ttdel)
            gamble_trial(ttdel) = [];
            Rtime(ttdel) = [];
        end
               
    
        if isempty(data_outcome_all)
%             ntrials_outcome{sn} = 1:length(data_outcome.time);
%             ntrials_pRPE{sn} = 1:length(data_positive.time);
%             ntrials_nRPE{sn} = 1:length(data_negative.time);
            
            gamble_trials = find(gamble_trial);
            gamble_choices = find(gamble_choice);
            reaction_time = Rtime;
%             data_button_all = data_button;
            data_outcome_all = data_outcome;
            data_trial_all = data_trial;
            
            data_choice_all = data_choice;
            RPE_out = RPE;
            outcome_out = outcome;
        else
%             ntrials_outcome{sn} = cat(2,ntrials_outcome{sn},size(data_outcome_all.time,2)+[1:length(data_outcome.time)]);
%             ntrials_pRPE{sn} = cat(2,ntrials_pRPE{sn},size(data_pRPE.time,2)+[1:length(data_positive.time)]);
%             ntrials_nRPE{sn} = cat(2,ntrials_nRPE{sn},size(data_nRPE.time,2)+[1:length(data_negative.time)]);
       
%             data_button_all.time((end+1):(end+length(data_button.time))) = data_button.time;
%             data_button_all.trial((end+1):(end+length(data_button.time))) = data_button.trial;
%             data_button_all.sampleinfo((end+1):(end+length(data_button.time)),:) = ...
%                 data_button.sampleinfo;
            gamble_trials = cat(2,gamble_trials,find(gamble_trial)+size(gamble_trials,2));
            gamble_choices = cat(2,gamble_choices,find(gamble_choice)+size(gamble_choices,2));
            reaction_time = cat(2,reaction_time, Rtime);
            
            data_outcome_all.time((end+1):(end+length(data_outcome.time))) = data_outcome.time;
            data_outcome_all.trial((end+1):(end+length(data_outcome.time))) = data_outcome.trial;
            data_outcome_all.sampleinfo((end+1):(end+length(data_outcome.time)),:) = ...
                data_outcome.sampleinfo;
            RPE_out((end+1):(end+length(data_outcome.time))) = RPE;
            outcome_out((end+1):(end+length(data_outcome.time))) = outcome;
            
            data_trial_all.time((end+1):(end+length(data_trial.time))) = data_trial.time;
            data_trial_all.trial((end+1):(end+length(data_trial.time))) = data_trial.trial;
            data_trial_all.sampleinfo((end+1):(end+length(data_trial.time)),:) = ...
                data_trial.sampleinfo;
            
            data_choice_all.time((end+1):(end+length(data_choice.time))) = data_choice.time;
            data_choice_all.trial((end+1):(end+length(data_choice.time))) = data_choice.trial;
            data_choice_all.sampleinfo((end+1):(end+length(data_choice.time)),:) = ...
                data_choice.sampleinfo;
            
           
        end
        
        
    end
      
    cfg = [];
    cfg.trials = gamble_trials;
    timelock_trials{sn} = ft_timelockanalysis(cfg,data_trial_all); 
    RT_gamble{sn} = reaction_time(gamble_trials);   
    
    cfg.trials = 1:size(data_trial_all.time,2); cfg.trials(gamble_trials) = [];
    timelock_trialsc{sn} = ft_timelockanalysis(cfg,data_trial_all); 
    RT_certain{sn} = reaction_time(cfg.trials);
    
    cfg = [];
    cfg.trials = gamble_choices;
    timelock_choice{sn} = ft_timelockanalysis(cfg,data_choice_all);
 
    cfg = [];
    cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
    timelock_choicec{sn} = ft_timelockanalysis(cfg,data_choice_all);
    
    timelock_out{sn} = ft_timelockanalysis([],data_outcome_all);
    
    [RPE,ind] = sort(RPE_out);
    outcome_sort = outcome_out(ind);
    pind = ind(RPE>5);
    nind = ind(RPE<-5 & outcome_sort< 0);
    npos = length(pind);
    nneg = length(nind);
    
    if npos > nneg
        pind = pind((end-nneg+1):end);
    elseif npos < nneg 
        keyboard
        [outcome_sort,ind] = sort(outcome_out,'descend');
        RPE = RPE_out(ind);
        pind = ind(RPE>1.5);
        pind = pind(1:nneg);
    end
    
    cfg = [];
    cfg.trials = pind;
    timelock_pRPE{sn} = ft_timelockanalysis(cfg,data_outcome_all);
    cfg.trials = nind;
    timelock_nRPE{sn} = ft_timelockanalysis(cfg,data_outcome_all);
    
    % TFS
    cfg  =[];
    cfg.taper = 'hanning';
    cfg.output = 'pow';
    cfg.pad = 'nextpow2';
    cfg.method = 'mtmconvol';
%     cfg.foi    = 2:2:34;
%     cfg.t_ftimwin = ones(length(cfg.foi),1)*0.5;
    
    cfg.foi          = 2:1:30;
    cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window    
    
    cfg.keeptrials = 'no';
    
    cfg.toi    = outcome_wind(1):0.05:outcome_wind(2);
    tfs_out{sn}    = ft_freqanalysis(cfg,data_outcome_all);
    
    cfg.trials = pind;
    tfs_pRPE{sn}   = ft_freqanalysis(cfg,data_outcome_all);
    cfg.trials = nind;
    tfs_nRPE{sn}   = ft_freqanalysis(cfg,data_outcome_all);
    
    cfg.trials  = gamble_choices;
    cfg.toi = choice_wind(1):0.05:choice_wind(2);
    tfs_choice{sn}   = ft_freqanalysis(cfg,data_choice_all);
    
    cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
    tfs_choicec{sn}   = ft_freqanalysis(cfg,data_choice_all);
    
    cfg.toi = outcome_wind(1):0.05:outcome_wind(2);
    cfg.trials  = gamble_trials;
    tfs_trials{sn}   = ft_freqanalysis(cfg,data_trial_all);
    
    cfg.trials = 1:size(data_trial_all.time,2); cfg.trials(gamble_trials) = [];
    tfs_trialsc{sn} = ft_freqanalysis(cfg,data_trial_all); 
    
end


%%
sn = [1:9,11,12,14:16];
% sn = [1:6,8,9,14];

grandavg_out = ft_timelockgrandaverage([],timelock_out{sn});
grandavg_pRPE = ft_timelockgrandaverage([],timelock_pRPE{sn});
grandavg_nRPE = ft_timelockgrandaverage([],timelock_nRPE{sn});
grandavg_trials = ft_timelockgrandaverage([],timelock_trials{sn});
grandavg_trialsc = ft_timelockgrandaverage([],timelock_trialsc{sn});
grandavg_choice = ft_timelockgrandaverage([],timelock_choice{sn});
grandavg_choicec = ft_timelockgrandaverage([],timelock_choicec{sn});

tfsavg_out = ft_freqgrandaverage([],tfs_out{sn});
tfsavg_pRPE = ft_freqgrandaverage([],tfs_pRPE{sn});
tfsavg_nRPE = ft_freqgrandaverage([],tfs_nRPE{sn});
tfsavg_trials = ft_freqgrandaverage([],tfs_trials{sn});
tfsavg_trialsc = ft_freqgrandaverage([],tfs_trialsc{sn});
tfsavg_choice = ft_freqgrandaverage([],tfs_choice{sn});
tfsavg_choicec = ft_freqgrandaverage([],tfs_choicec{sn});


%% All outcomes big peaks
% 

nsens = 5;
times = [0.10 0.12; 0.23 0.3; 0.38 0.42; 0.55 0.61];
% find sensors with highest ER
figure(1); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_out.time>= timee(1) & grandavg_out.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_out.avg(:,timesind) < (-3e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
[~,ind] = sort(mean(grandavg_out.avg(:,timesind),2 ));
sensors = data_outcome_all.label(ind(1:nsens));
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1 1]*0.5e-13;
cfg.channel = sensors;
subplot(431)
ft_singleplotER(cfg,grandavg_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
ylabel('ERF (T)')

subplot(432)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (3e-14),2 ));
cfgf = [];
cfgf.baseline     = [-1 -.05];
cfgf.baselinetype = 'absolute';
cfgf.maskstyle   = 'saturation';
cfgf.interactive  = 'no';
cfgf.zlim         = [-1 1]*5e-28;
cfgf.xlim         = [-1 3];
cfgf.ylim        = [1 35];
cfgf.trials  = 'all';
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_out);
ylabel('Frequency (Hz)')

% sensors = 'MRT*';
timee = times(2,:);
timesind = grandavg_out.time>= timee(1) & grandavg_out.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_out.avg(:,timesind) < (-5e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_out.avg(:,timesind),2 ));
sensors = data_outcome_all.label(ind(1:nsens));
cfg.channel = sensors;
subplot(434)
ft_singleplotER(cfg,grandavg_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
ylabel('ERF (T)')

subplot(435)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (5e-14),2 ));
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_out);
ylabel('Frequency (Hz)')


% sensors = 'MRC*';
timee = times(3,:);
timesind = grandavg_out.time>= timee(1) & grandavg_out.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_out.avg(:,timesind) > (4e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_out.avg(:,timesind),2 ),'descend');
sensors = data_outcome_all.label(ind(1:nsens));
cfg.channel = sensors;
subplot(437)
ft_singleplotER(cfg,grandavg_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
ylabel('ERF (T)')

subplot(438)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (5e-14),2 ));
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_out);
ylabel('Frequency (Hz)')

% sensors = 'MRP*';
timee = times(4,:);
timesind = grandavg_out.time>= timee(1) & grandavg_out.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_out.avg(:,timesind) > (3e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_out.avg(:,timesind),2 ),'descend');
sensors = data_outcome_all.label(ind(1:nsens));
cfg.channel = sensors;
subplot(4,3,10)
ft_singleplotER(cfg,grandavg_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
xlabel('Time (s)'); ylabel('ERF (T)')


subplot(4,3,11)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (4e-14),2 ));
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_out);
ylabel('Frequency (Hz)'); xlabel('Time (s)')


for ii = 1:size(times,1)
    subplot(4,3,(ii-1)*3+3)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,grandavg_out);
end
set(gcf,'color','w',...
    'name',sprintf('Gamble outcome at t = 0 (%d trials)',length(data_outcome_all.time)))

%% Plot TFS on sensors
time1 = -0.5:0.5:2;
time2 = 0:0.5:2.5;

freq = [4 8; 8 13; 14 23];
z = [4 6 2]*1e-28;
for jj = 1:3
    figure(jj+1); clf
    for ii = 1:length(time1)
        subplot(2,3,ii)
        cfg = [];
        cfg.baseline     = [-1 -0.05];
        cfg.baselinetype = 'absolute';
        cfg.xlim         = [time1(ii) time2(ii)];
        cfg.zlim         = [-1 1]*z(jj);
        cfg.ylim         = freq(jj,:);
        cfg.marker       = 'on';
        cfg.layout = 'CTF275_helmet.mat';
        ft_topoplotTFR(cfg, tfs_pRPE{16});
    end
end

%% Positive and negative outcomes

tfsbase_out = zeros(size(tfsavg_out.label,1),size(tfsavg_out.freq,2));
for ii = 1:size(tfsavg_out.freq,2)
    ind = find(squeeze(~isnan(tfsavg_out.powspctrm(1,ii,1:20))));
    if ~isempty(ind)
        tfs = tfsavg_out.powspctrm(:,ii,ind(1):20);    
        tfsbase_out(:,ii) =  mean(tfs,3);
    end
end
grandavg_diff = grandavg_pRPE;
grandavg_diff.avg = grandavg_nRPE.avg - grandavg_pRPE.avg;
times = [0.2 0.3];

% find sensors with highest ER
figure(3); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_out.time>= timee(1) & grandavg_out.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_out.avg(:,timesind) > (6e-14),2 ));
% sensors = sensors(strncmp(sensors,'MR',2)); % limit to right hemisphere
[~,ind] = sort(mean(grandavg_out.avg(:,timesind),2 ),'descend');
sensors = data_outcome_all.label(ind(1:nsens));
sensorstfs = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (6e-14),2 ));
sensorstfs = sensors;
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(431)
ft_singleplotER(cfg,grandavg_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('All outcomes grand average')

cfg.channel = sensors;
subplot(434)
ft_singleplotER(cfg,grandavg_pRPE);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Positive outcome')

subplot(437)
ft_singleplotER(cfg,grandavg_nRPE);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative outcome')

cfg.channel = sensors;
subplot(4,3,10)
ft_singleplotER(cfg,grandavg_diff);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative - Positive outcome')
xlabel('Time (s)')
  
cfgf = [];
cfgf.baseline  = [-2 -0.1];
cfgf.baselinetype = 'absolute';
freq = ft_freqbaseline(cfgf, tfsavg_out);

cfgf = [];
% cfgf.baseline     = [-2 -0.1];
% cfgf.baselinetype = 'absolute';
cfgf.maskstyle   = 'saturation';
cfgf.interactive  = 'no';
cfgf.zlim         = [-1 1]*0.5e-27;
cfgf.xlim         = [-0.7 2.7];
cfgf.ylim        = [1 35];
cfgf.trials  = 'all';
cfgf.channel  =  sensorstfs;

subplot(432)
tfs_plot = tfsavg_out;
tfs_plot.powspctrm = tfs_plot.powspctrm - tfsbase_out;
ft_singleplotTFR(cfgf, tfs_plot);

subplot(435)
tfs_plot = tfsavg_pRPE;
tfs_plot.powspctrm = tfs_plot.powspctrm - tfsbase_out;
ft_singleplotTFR(cfgf, tfs_plot);

subplot(438)
tfs_plot = tfsavg_nRPE;
tfs_plot.powspctrm = tfs_plot.powspctrm - tfsbase_out;
ft_singleplotTFR(cfgf, tfs_plot);

tfs_diff = tfsavg_nRPE;
tfs_diff.powspctrm = tfsavg_nRPE.powspctrm - tfsavg_pRPE.powspctrm;
subplot(4,3,11)
ft_singleplotTFR(cfgf, tfs_diff);



for ii = 1:4
    subplot(4,3,(ii-1)*3+3)
    cfg = [];
    cfg.xlim = times;
    cfg.zlim = [-1 1]*5e-14;
    cfg.layout = 'CTF275_helmet.mat';
    switch ii
        case 1
            ft_topoplotTFR(cfg,grandavg_out);
        case 2
            ft_topoplotTFR(cfg,grandavg_pRPE);
        case 3
            ft_topoplotTFR(cfg,grandavg_nRPE);
        case 4
            ft_topoplotTFR(cfg,grandavg_diff);
    end
end

set(gcf,'color','w',...
    'name',sprintf('Gamble outcome at t = 0'))


%% All choices (after decision to gamble or not)
fname = sprintf('Chosen gamble option at t=0 ');
times = [-0.1 0.02 ; 0.02 0.13; 0.13 0.2; 0.25 0.35];


% find sensors with highest ER
figure(4); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_choice.time>= timee(1) & grandavg_choice.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_choice.avg(:,timesind) < (-6e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
[~,ind] = sort(mean(grandavg_choice.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-1 1.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(431)
ft_singleplotER(cfg,grandavg_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

subplot(432)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (3e-14),2 ));
cfgf = [];
cfgf.baseline     = [-2 -1];
cfgf.baselinetype = 'absolute';
cfgf.maskstyle   = 'saturation';
cfgf.interactive  = 'no';
cfgf.zlim         = [-1 1]*5e-28;
cfgf.xlim         = [-2.5 2.5];
cfgf.ylim        = [1 35];
cfgf.trials  = 'all';
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_choice);
ylabel('Frequency (Hz)')


% sensors = 'MRT*';
timee = times(2,:);
timesind = grandavg_choice.time>= timee(1) & grandavg_choice.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_choice.avg(:,timesind) > (5e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_choice.avg(:,timesind),2 ),'descend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
subplot(434)
ft_singleplotER(cfg,grandavg_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

subplot(435)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_choice);
ylabel('Frequency (Hz)')


% sensors = 'MRC*';
timee = times(3,:);
timesind = grandavg_choice.time>= timee(1) & grandavg_choice.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_choice.avg(:,timesind) < (-4e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_choice.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
subplot(437)
ft_singleplotER(cfg,grandavg_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

subplot(438)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_choice);
ylabel('Frequency (Hz)')

% sensors = 'MRP*';
timee = times(4,:);
timesind = grandavg_choice.time>= timee(1) & grandavg_choice.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_choice.avg(:,timesind) < (-5e-14),2 ));
% sensors = sensors(strncmp(sensors,'MR',2));
[~,ind] = sort(mean(grandavg_choice.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
subplot(4,3,10)
ft_singleplotER(cfg,grandavg_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
xlabel('Time (s)')

subplot(4,3,11)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_choice);
ylabel('Frequency (Hz)')


for ii = 1:size(times,1)
    subplot(4,3,(ii-1)*3+3)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,grandavg_choice);
end
set(gcf,'color','w','name',fname);

%% All choices (after decision not to gamble)
cfg = [];
cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
% grandavg_choicec = ft_grandavganalysis(cfg,data_choice_all);
fname = sprintf('Chosen certain option at t=0 (%d trials)',length(cfg.trials));
% find sensors with highest ER
figure(5); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_choicec.time>= timee(1) & grandavg_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_choicec.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-1 1.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,grandavg_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRT*';
timee = times(2,:);
timesind = grandavg_choicec.time>= timee(1) & grandavg_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_choicec.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,grandavg_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRC*';
timee = times(3,:);
timesind = grandavg_choicec.time>= timee(1) & grandavg_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_choicec.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,grandavg_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRP*';
timee = times(4,:);
timesind = grandavg_choicec.time>= timee(1) & grandavg_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_choicec.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,grandavg_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
xlabel('Time (s)')

for ii = 1:size(times,1)
    subplot(4,2,(ii-1)*2+2)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,grandavg_choicec);
end

cfg = [];
cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
set(gcf,'color','w','name',fname);

%% All trials (option presentation)

fname = sprintf('Options at t=0, chosen gamble at t in red');

b = 15; % Number of bins for reaction time distribution

RT = cell2mat(RT_gamble);

times = [ 0.07 0.15; 0.18 0.25; 0.3 0.4;0.45 0.8];
% find sensors with highest ER
figure(6); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_trials.time>= timee(1) & grandavg_trials.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_trials.avg(:,timesind) < (-5e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
[~,ind] = sort(mean(grandavg_trials.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));


% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 2.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(431)
ft_singleplotER(cfg,grandavg_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

subplot(432)
% sensors = data_outcome_all.label(any(abs(grandavg_out.avg(:,timesind)) > (3e-14),2 ));
cfgf = [];
cfgf.baseline     = [-1 -0.05];
cfgf.baselinetype = 'absolute';
cfgf.maskstyle   = 'saturation';
cfgf.interactive  = 'no';
cfgf.zlim         = [-1 1]*5e-28;
cfgf.xlim         = [-1 3];
cfgf.ylim        = [1 35];
cfgf.trials  = 'all';
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_trials);
ylabel('Frequency (Hz)')


% sensors = 'MRT*';
timee = times(2,:);
timesind = grandavg_trials.time>= timee(1) & grandavg_trials.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_trials.avg(:,timesind) > (4e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_trials.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
subplot(434)
ft_singleplotER(cfg,grandavg_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

subplot(435)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_trials);
ylabel('Frequency (Hz)')


% sensors = 'MRC*';
timee = times(3,:);
timesind = grandavg_trials.time>= timee(1) & grandavg_trials.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_trials.avg(:,timesind) < (-4e-14),2 ));
% sensors = sensors(strncmp(sensors,'MR',2));
[~,ind] = sort(mean(grandavg_trials.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
subplot(437)
ft_singleplotER(cfg,grandavg_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

subplot(438)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_trials);
ylabel('Frequency (Hz)')

% sensors = 'MRP*';
timee = times(4,:);
timesind = grandavg_trials.time>= timee(1) & grandavg_trials.time<=timee(2);
% sensors = data_outcome_all.label(any(grandavg_trials.avg(:,timesind) < (-6e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
[~,ind] = sort(mean(grandavg_trials.avg(:,timesind),2 ),'ascend');
sensors = data_outcome_all.label(ind(1:nsens));

cfg.channel = sensors;
cfg.ylim = [-1 1]*2e-13;
subplot(4,3,10)
ft_singleplotER(cfg,grandavg_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

xlabel('Time (s)')

subplot(4,3,11)
cfgf.channel  =  sensors;
ft_singleplotTFR(cfgf, tfsavg_trials);
ylabel('Frequency (Hz)');xlabel('Time (s)')



for ii = 1:size(times,1)
    subplot(4,3,(ii-1)*3+3)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,grandavg_trials);
end
set(gcf,'color','w','name',fname);

%%
cfg= [];
cfg.trials = 1:size(data_trial_all.time,2); cfg.trials(gamble_trials) = [];
fname = sprintf('Options at t=0, chosen certain at t in red (%d trials)',length(cfg.trials));

RT = reaction_time;
RT(gamble_trials) = [];

% grandavg_trialsc = ft_grandavganalysis(cfg,data_trial_all);
% find sensors with highest ER
figure(7); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = grandavg_trialsc.time>= timee(1) & grandavg_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_trialsc.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 2.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,grandavg_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRT*';
timee = times(2,:);
timesind = grandavg_trialsc.time>= timee(1) & grandavg_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_trialsc.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,grandavg_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRC*';
timee = times(3,:);
timesind = grandavg_trialsc.time>= timee(1) & grandavg_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_trialsc.avg(:,timesind) < (-4e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,grandavg_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRP*';
timee = times(4,:);
timesind = grandavg_trialsc.time>= timee(1) & grandavg_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(grandavg_trialsc.avg(:,timesind) < (-6e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,grandavg_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

xlabel('Time (s)')

for ii = 1:size(times,1)
    subplot(4,2,(ii-1)*2+2)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,grandavg_trialsc);
end
set(gcf,'color','w','name',fname);
