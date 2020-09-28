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

data_button_all = [];
data_outcome_all = [];
data_trial_all = [];
data_choice_all = [];
data_pRPE = [];
data_nRPE = [];

ntrials_pRPE = cell(1,16);
ntrials_nRPE = cell(1,16);
ntrials_outcome = cell(1,16);
ntrials_button = cell(1,16);

data_lock_diff = cell(1,16);
data_lock_outcome = cell(1,16);
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
for sn = [1:6,8,9,14] % [1:9,11,12,14:16] 
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
            load([processing_folder,'/ICA_artifacts.mat']);
            cfg           = [];
            cfg.component = 1:length(comps.label);
            data_clean    = ft_rejectcomponent(cfg, comps,data);
        end
        
        data = data_clean;
        clear data_clean
        
        filt_order = []; % default
        % data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, freq_band(freq,:),filt_order,'but');
        data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 ); % 35Hz for ERF
        
        data.trial{1} = data_filt;
       
        clear data_filt
        
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
            
        [RPE,ind] = sort(outcome_match.RPE,'descend');
        pind = ind(RPE>4);
        [RPE,ind] = sort(outcome_match.RPE,'ascend');
        nind = ind(RPE<-4);
        
        data_positive = define_trials(outcome_match.sample(pind),data,bv_match.time,[-2 3]);
        data_negative = define_trials(outcome_match.sample(nind),data,bv_match.time,[-2 3]);
        
        npos = size(data_positive.trial,2);
        nneg = size(data_negative.trial,2);
        if npos > nneg
            data_positive.trial = data_positive.trial(1:nneg);
            data_positive.time = data_positive.time(1:nneg);
            data_positive.sampleinfo = data_positive.sampleinfo(1:nneg,:);
        elseif npos < nneg
            data_negative.trial = data_negative.trial(1:npos);
            data_negative.time = data_negative.time(1:npos);
            data_negative.sampleinfo = data_negative.sampleinfo(1:npos,:);
        end

%         data_button = define_trials(buttonpress,data,bv_match.time,[-1 2]);
        
        Rsamples = outcome_match.sample(outcome_match.win~=0);
        data_outcome = define_trials(Rsamples,data,bv_match.time,[-2 3]);
        
        % Select trials where participant choose an option
        
        Rsamples = choice_match.sample(choice_match.sample~=0 & ~isnan(choice_match.choice));
        [data_choice,ttdel] = define_trials(Rsamples,data,bv_match.time,[-2 2]);
        
        gamble_choice =  choice_match.choice(choice_match.sample~=0 & ~isnan(choice_match.choice)) == 1;
        if ~isempty(ttdel)
            gamble_choice(ttdel) = [];
        end
        
        Rsamples = answer_match.sample(answer_match.sample~=0 & ~isnan(answer_match.choice));
        [data_trial,ttdel] = define_trials(Rsamples,data,bv_match.time,[-1 2]);
        
        gamble_trial = answer_match.choice(answer_match.sample~=0 & ~isnan(answer_match.choice)) == 1;
        Rtime = answer_match.RT(answer_match.sample~=0 & ~isnan(answer_match.choice));
        if ~isempty(ttdel)
            gamble_trial(ttdel) = [];
            Rtime(ttdel) = [];
        end
       
        data_outlock = ft_timelockanalysis([],data_outcome);
        clear vol
        load([processing_folder,'/headmodel.mat']);
        cfg = [];
        cfg.inwardshift  = 10; %mm
        cfg.headmodel = vol;
        cfg.method  = 'template';
        % transform data to planar gradiometer
        cfg.feedback        = 'no';
        cfg.neighbours      = ft_prepare_neighbours(cfg, data_outlock);
        cfg.planarmethod    = 'sincos';
              
        avgFICplanar        = ft_megplanar(cfg, data_outlock);
        avgFICplanarComb   = ft_combineplanar([],avgFICplanar);
        avgFICplanarComb.dof = data_outlock.dof;
    
        data_lock_outcome{sn} = avgFICplanarComb;
        
        timelock_pRPE = ft_timelockanalysis([],data_positive);
        timelock_nRPE = ft_timelockanalysis([],data_negative);

        timelock_diff = timelock_nRPE;
        timelock_diff.avg = timelock_nRPE.avg - timelock_pRPE.avg;
        
        avgFICplanar        = ft_megplanar(cfg, timelock_diff);
        avgFICplanarComb   = ft_combineplanar([],avgFICplanar);
        avgFICplanarComb.dof = timelock_diff.dof;
    
        data_lock_diff{sn} = avgFICplanarComb;
        
        
        if isempty(data_outcome_all)
            ntrials_outcome{sn} = 1:length(data_outcome.time);
            ntrials_pRPE{sn} = 1:length(data_positive.time);
            ntrials_nRPE{sn} = 1:length(data_negative.time);
            
%             ntrials_button{sn}  = 1:length(data_button.time);
            
            gamble_trials = find(gamble_trial);
            gamble_choices = find(gamble_choice);
            reaction_time = Rtime;
%             data_button_all = data_button;
            data_outcome_all = data_outcome;
            data_trial_all = data_trial;
            data_pRPE = data_positive;
            data_nRPE = data_negative;
            data_choice_all = data_choice;
        else
            ntrials_outcome{sn} = cat(2,ntrials_outcome{sn},size(data_outcome_all.time,2)+[1:length(data_outcome.time)]);
            ntrials_pRPE{sn} = cat(2,ntrials_pRPE{sn},size(data_pRPE.time,2)+[1:length(data_positive.time)]);
            ntrials_nRPE{sn} = cat(2,ntrials_nRPE{sn},size(data_nRPE.time,2)+[1:length(data_negative.time)]);
%             ntrials_button{sn}  = cat(2,ntrials_button{sn},size(data_button.time,2)+[1:length(data_button.time)]);
            
%             data_button_all.time((end+1):(end+length(data_button.time))) = data_button.time;
%             data_button_all.trial((end+1):(end+length(data_button.time))) = data_button.trial;
%             data_button_all.sampleinfo((end+1):(end+length(data_button.time)),:) = ...
%                 data_button.sampleinfo;
            gamble_trials = cat(2,gamble_trials,find(gamble_trial)+size(gamble_trials,2));
            gamble_choices = cat(2,gamble_choices,find(gamble_choice)+size(gamble_choices,2));
            reaction_time = cat(2,reaction_time, Rtime);
%             
            data_outcome_all.time((end+1):(end+length(data_outcome.time))) = data_outcome.time;
            data_outcome_all.trial((end+1):(end+length(data_outcome.time))) = data_outcome.trial;
            data_outcome_all.sampleinfo((end+1):(end+length(data_outcome.time)),:) = ...
                data_outcome.sampleinfo;
            
            data_trial_all.time((end+1):(end+length(data_trial.time))) = data_trial.time;
            data_trial_all.trial((end+1):(end+length(data_trial.time))) = data_trial.trial;
            data_trial_all.sampleinfo((end+1):(end+length(data_trial.time)),:) = ...
                data_trial.sampleinfo;
            
            data_choice_all.time((end+1):(end+length(data_choice.time))) = data_choice.time;
            data_choice_all.trial((end+1):(end+length(data_choice.time))) = data_choice.trial;
            data_choice_all.sampleinfo((end+1):(end+length(data_choice.time)),:) = ...
                data_choice.sampleinfo;
            
            
            data_pRPE.time((end+1):(end+length(data_positive.time))) = data_positive.time;
            data_pRPE.trial((end+1):(end+length(data_positive.time))) = data_positive.trial;
            data_pRPE.sampleinfo((end+1):(end+length(data_positive.time)),:) = ...
                data_positive.sampleinfo;
            
            
            data_nRPE.time((end+1):(end+length(data_negative.time))) = data_negative.time;
            data_nRPE.trial((end+1):(end+length(data_negative.time))) = data_negative.trial;
            data_nRPE.sampleinfo((end+1):(end+length(data_negative.time)),:) = ...
                data_negative.sampleinfo;
%             
        end
        
        
    end
    
    %         save([data_path,data_names{runs},'/beamforming/VE_erf.mat'],'VEp','VEn','Rwin')
end

%%

timelock_pRPE = ft_timelockanalysis([],data_pRPE);
timelock_nRPE = ft_timelockanalysis([],data_nRPE);

timelock_diff = timelock_nRPE;
timelock_diff.avg = timelock_nRPE.avg - timelock_pRPE.avg;
% time lock average over all trials and subjects in axial gradiometers
timelock_out = ft_timelockanalysis([],data_outcome_all);


% Grand average of planar gradiometers
cfg   =[];
cfg.method = 'across'; %'within'
[grandavg_outcome] = ft_timelockgrandaverage(cfg, data_lock_outcome{1}, data_lock_outcome{2}, data_lock_outcome{3},...
    data_lock_outcome{4}, data_lock_outcome{5}, data_lock_outcome{6}, data_lock_outcome{8}, data_lock_outcome{9}, data_lock_outcome{14});

cfg  =[];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.toi    = -2:0.05:3;
cfg.foi    = 1:1:35;
cfg.width = 7;
cfg.keeptrials = 'no';
outcome_tfs    = ft_freqanalysis(cfg,data_outcome_all);


% Grand average of planar gradiometers
cfg   =[];
cfg.method = 'across'; %'within'
[grandavg_diff] = ft_timelockgrandaverage(cfg, data_lock_diff{1}, data_lock_diff{2}, data_lock_diff{3},...
    data_lock_diff{4}, data_lock_diff{5}, data_lock_diff{6}, data_lock_diff{8}, data_lock_diff{9}, data_lock_diff{14});

cfg  =[];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.toi    = -2:0.05:3;
cfg.foi    = 1:1:35;
cfg.width = 3;
cfg.keeptrials = 'no';
pRPE_tfs    = ft_freqanalysis(cfg,data_pRPE);
nRPE_tfs    = ft_freqanalysis(cfg,data_nRPE);

diff_tfs   = nRPE_tfs;
diff_tfs.powspctrm = nRPE_tfs.powspctrm  - pRPE_tfs.powspctrm;
%% Plots of outcome
% times = [0.03 0.04; 0.125 0.135;]; % button presses
times = [0.22 0.3; 0.23 0.28; 0.37 0.42; 0.4 0.6];

for ii =1:size(times,1)
    cfg = [];
    figure(ii)
    clf
    subplot(421);
    cfg.xlim = times(ii,:);
    cfg.zlim = [-7 7]*1e-14;
    cfg.colorbar = 'no';
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotER(cfg,timelock_out)
    
    subplot(322);
    cfg.zlim = [-7 7]*1e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotER(cfg,grandavg_outcome);
    colormap jet
    colorbar
    subplot(4,2,3:4)
    
    timee = times(ii,:);
    timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);   
    
%     sensors = data_outcome_all.label(any( (grandavg.avg(:,timesind)) > (3e-14),2 )); 
    [~,ind] = sort(mean(abs(grandavg_outcome.avg(:,timesind)) ,2 ) ,'descend') ;
    sensors = data_outcome_all.label(ind(1:7));
    if ii == 1
        sensors = data_outcome_all.label(any( (grandavg_outcome.avg(:,timesind)) > (4.2e-14),2 )); 
        sensors = sensors(strncmp(sensors,'MRO',3) | strncmp(sensors,'MLO',3));
    end
  
    cfg = [];
    cfg.xlim = [-1 2];
    cfg.ylim = [0 1]*1e-13;
    cfg.channel = sensors;
    ft_singleplotER(cfg,grandavg_outcome);
    hold on 
    fill([timee(1) timee(2) timee(2) timee(1)],[0 0 1 1]*1e-13,[1 0 0],'facealpha',0.2,'edgecolor','none')
   
      
    subplot(4,2,5:6)
    
    timee = times(ii,:);
    
    timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
%     sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) < (-5e-14),2 ));
%     sensors = sensors(strncmp(sensors,'ML',2));

    [~,ind] = sort(mean((timelock_out.avg(:,timesind)) ,2 ) ,'ascend') ;
    sensors = data_outcome_all.label(ind(1:7));
    if ii == 1
        sensors = data_outcome_all.label(any( (timelock_out.avg(:,timesind)) > (4e-14),2 )); 
        sensors = sensors(strncmp(sensors,'MLO',3) );
    end
  
    cfg = [];
    cfg.xlim = [-1 2];
    cfg.ylim = [-0.8 0.6]*1e-13;
    cfg.channel = sensors;
    ft_singleplotER(cfg,timelock_out);
    hold on 
    fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1]*1e-13,[1 0 0],'facealpha',0.2,'edgecolor','none')
   
    subplot(4,2,7:8)
    timee = times(ii,:);
    timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
%     sensors = data_outcome_all.label(any( abs(timelock_outcome.avg(:,timesind)) > (6e-14),2 )); 
    [~,ind] = sort(mean(abs(timelock_out.avg(:,timesind)) ,2 ) ,'descend') ;
    sensors = data_outcome_all.label(ind(1:7));
    if ii == 1
        sensors = data_outcome_all.label(any( abs(timelock_out.avg(:,timesind)) > (4e-14),2 )); 
        sensors = sensors(strncmp(sensors,'MRO',3) | strncmp(sensors,'MLO',3));
    end
    
    cfg = [];
    cfg.baseline     = [-1 -0.1];
    cfg.baselinetype = 'absolute';
    cfg.maskstyle   = 'saturation';
    cfg.interactive  = 'no';
    cfg.zlim         = [-1 1]*8e-25;
    cfg.xlim         = [-1 2];
    cfg.ylim        = [2 25];
    cfg.trials  = 'all';
    cfg.channel  =  sensors;
    ft_singleplotTFR(cfg, outcome_tfs);
    colorbar off
    set(gcf,'color','w','name',sprintf('Time %.2f-%.2f',timee(1),timee(2)))
end

%% Plots of positive and negative RPE
% times = [0.03 0.04; 0.125 0.135;]; % button presses

% cfg = [];
% cfg.inwardshift  = 10; %mm
% % cfg.headmodel = vol;
% cfg.method  = 'template';
% % transform data to planar gradiometer
% cfg.feedback        = 'no';
% cfg.neighbours      = ft_prepare_neighbours(cfg, timelock_pRPE);
% cfg.planarmethod    = 'sincos';
% avgFICplanar        = ft_megplanar(cfg, timelock_pRPE);
% grandavg_pRPE   = ft_combineplanar([],avgFICplanar);
% 
% avgFICplanar        = ft_megplanar(cfg, timelock_nRPE);
% grandavg_nRPE   = ft_combineplanar([],avgFICplanar);

times = [0.18 0.3; 0.22 0.3; 0.25 0.3; 0.37 0.42;];

ii =1;
timee = times(ii,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);

[~,ind] = sort(mean(abs(grandavg_outcome.avg(:,timesind)) ,2 ) ,'descend') ;
sensorsp = data_outcome_all.label(ind(1:7));

[~,ind] = sort(mean((timelock_out.avg(:,timesind)) ,2 ) ,'ascend') ;
sensors = data_outcome_all.label(ind(1:7));

for jj = 1:3
    switch jj
        case 1
            timelock_data = timelock_pRPE;
            data_tfs = pRPE_tfs;
            grandavg = grandavg_pRPE;
        case 2 
            timelock_data = timelock_nRPE;
            data_tfs = nRPE_tfs;
            grandavg = grandavg_nRPE;
        case 3
            timelock_data = timelock_diff;
            data_tfs = diff_tfs;  
            grandavg = grandavg_diff;
            grandavg.avg = grandavg_nRPE.avg - grandavg_pRPE.avg;
    end
    
    ii = 2;
    cfg = [];
    figure(jj+4)
    clf
    subplot(421);
    cfg.xlim = times(ii,:);
    cfg.zlim = [-7 7]*1e-14;
    if jj == 3
        cfg.zlim = [-5 5]*1e-14;
    end
    cfg.colorbar = 'no';
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotER(cfg,timelock_data)
    
    subplot(422);
    cfg.zlim = [-8 8]*1e-14;
    if jj == 3
        cfg.zlim = [-4 4]*1e-14;
    end
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotER(cfg,grandavg);
    colormap jet
    
    subplot(4,2,3:4)
    
    timee = times(ii,:);
    timesind = timelock_pRPE.time>= timee(1) & timelock_pRPE.time<=timee(2);   
    
%     [~,ind] = sort(mean(abs(grandavg_outcome.avg(:,timesind)) ,2 ) ,'descend') ;
%     sensors = data_outcome_all.label(ind(1:5));
%     if ii == 1
%         sensors = data_outcome_all.label(any( (grandavg_outcome.avg(:,timesind)) > (4.2e-14),2 )); 
%         sensors = sensors(strncmp(sensors,'MRO',3) | strncmp(sensors,'MLO',3));
%     end
  
    cfg = [];
    cfg.xlim = [-0.5 1.3];
    cfg.ylim = [0 1]*1e-13;
    if jj == 3
        cfg.ylim = [-0.5 0.5]*1e-13;
    end
    cfg.channel = sensorsp;
    ft_singleplotER(cfg,grandavg);
    hold on 
    fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1]*1e-13,[1 0 0],'facealpha',0.2,'edgecolor','none')
   
    
    subplot(4,2,5:6)    
    timee = times(ii,:);
    
    timesind = timelock_pRPE.time>= timee(1) & timelock_pRPE.time<=timee(2);
  
    cfg = [];
    cfg.xlim = [-1 2];
    cfg.ylim = [-0.8 0.6]*1e-13;
    cfg.channel = sensors;
    ft_singleplotER(cfg,timelock_data);
    hold on 
    fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1]*1e-13,[1 0 0],'facealpha',0.2,'edgecolor','none')
    
    
    
    subplot(4,2,7:8)
    timee = times(ii,:);
    timesind = timelock_pRPE.time>= timee(1) & timelock_pRPE.time<=timee(2);
%     [~,ind] = sort(mean(abs(timelock_data.avg(:,timesind)) ,2 ) ,'descend') ;
%     sensors = data_outcome_all.label(ind(1:5));
%     if ii == 1
%         sensors = data_outcome_all.label(any( abs(timelock_data.avg(:,timesind)) > (4e-14),2 )); 
%         sensors = sensors(strncmp(sensors,'MRO',3) | strncmp(sensors,'MLO',3));
%     end
    
    cfg = [];
    cfg.baseline     = [-1 -0.1];
    cfg.baselinetype = 'absolute';
    cfg.maskstyle   = 'saturation';
    cfg.interactive  = 'no';
    cfg.zlim         = [-1 1]*8e-25;
    cfg.xlim         = [-0.5 1.3];
    cfg.ylim        = [2 25];
    cfg.trials  = 'all';
    cfg.channel  =  sensors;
    ft_singleplotTFR(cfg, data_tfs);
    colorbar off
    set(gcf,'color','w')
end

%% All outcomes big peaks

times = [0.10 0.12; 0.23 0.3; 0.38 0.42; 0.55 0.61];
% find sensors with highest ER
figure(1); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) < (-3e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1 1]*0.5e-13;
cfg.channel = sensors;
subplot(431)
ft_singleplotER(cfg,timelock_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
% 
% subplot(432)
% sens = zeros(size(sensors));
% for ii = 1:length(sensors)
%    sens(ii) = find(strcmp(sensors{ii},data_outcome_all.label )) ;
% end
% pcolor(linspace(-.2,1,240),(highpass+lowpass)/2, ...
%     mean(outcome_filt_sens(:,:,sens)-tfsbase(:,:,sens),3)); shading interp
% colorbar; caxis([-1 1]*4e-15)

% sensors = 'MRT*';
timee = times(2,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) < (-4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(434)
ft_singleplotER(cfg,timelock_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% subplot(435)
% sens = zeros(size(sensors));
% for ii = 1:length(sensors)
%    sens(ii) = find(strcmp(sensors{ii},data_outcome_all.label )) ;
% end
% pcolor(linspace(-.2,1,240),(highpass+lowpass)/2, ...
%     mean(outcome_filt_sens(:,:,sens)-tfsbase(:,:,sens),3)); shading interp
% colorbar; caxis([-1 1]*4e-15)


% sensors = 'MRC*';
timee = times(3,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(437)
ft_singleplotER(cfg,timelock_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% subplot(438)
% sens = zeros(size(sensors));
% for ii = 1:length(sensors)
%    sens(ii) = find(strcmp(sensors{ii},data_outcome_all.label )) ;
% end
% pcolor(linspace(-.2,1,240),(highpass+lowpass)/2, ...
%     mean(outcome_filt_sens(:,:,sens)-tfsbase(:,:,sens),3)); shading interp
% colorbar; caxis([-1 1]*4e-15)

% sensors = 'MRP*';
timee = times(4,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) > (3e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(4,3,10)
ft_singleplotER(cfg,timelock_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
xlabel('Time (s)')

% subplot(4,3,11)
% sens = zeros(size(sensors));
% for ii = 1:length(sensors)
%    sens(ii) = find(strcmp(sensors{ii},data_outcome_all.label )) ;
% end
% pcolor(linspace(-.2,1,240),(highpass+lowpass)/2, ...
%     mean(outcome_filt_sens(:,:,sens)-tfsbase(:,:,sens),3)); shading interp
% colorbar; caxis([-1 1]*4e-15)

for ii = 1:size(times,1)
    subplot(4,3,(ii-1)*3+3)
    cfg = [];
    cfg.xlim = times(ii,:);
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    ft_topoplotTFR(cfg,timelock_out);
end
set(gcf,'color','w',...
    'name',sprintf('Gamble outcome at t = 0 (%d trials)',length(data_outcome_all.time)))

%% Positive and negative outcomes

times = [0.18 0.3];
% find sensors with highest ER
figure(2); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1 1]*0.8e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_out);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('All outcomes grand average')

% sensors = 'MRT*';
% timee = times(2,:);
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) < (-5e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_pRPE);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Positive outcome')

% sensors = 'MRC*';
% timee = times(3,:);
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) > (4e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
% cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_nRPE);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative outcome')

% sensors = 'MRP*';
% timee = times(4,:);
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensors = data_outcome_all.label(any(timelock_out.avg(:,timesind) > (4e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_diff);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative - Positive outcome')
xlabel('Time (s)')

for ii = 1:4
    subplot(4,2,(ii-1)*2+2)
    cfg = [];
    cfg.xlim = times;
    cfg.zlim = [-1 1]*4e-14;
    cfg.layout = 'CTF275_helmet.mat';
    switch ii
        case 1
            ft_topoplotTFR(cfg,timelock_out);
        case 2
            ft_topoplotTFR(cfg,timelock_pRPE);
        case 3
            ft_topoplotTFR(cfg,timelock_nRPE);
        case 4
            ft_topoplotTFR(cfg,timelock_diff);
    end
end


%% One subject trials
sn =2;

cfg = [];
cfg.trials = ntrials_pRPE{sn};
timelock_positive = ft_timelockanalysis(cfg,data_pRPE);
cfg.trials = ntrials_nRPE{sn};
timelock_negative = ft_timelockanalysis(cfg,data_nRPE);

cfg.trials = ntrials_outcome{sn};
timelock_outcome = ft_timelockanalysis(cfg,data_outcome_all);

timelock_diff1 = timelock_positive;
timelock_diff1.avg = timelock_negative.avg - timelock_positive.avg;


times = [0.18 0.25];
% times = [0.25 0.35];
% times = [0.33 0.42];
% times = [0.3 0.4];
timee = times(1,:);
timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);

% sensors = data_outcome_all.label(any(timelock_outcome.avg(:,timesind) < (-5e-14),2 ));
% sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere

sensors = data_outcome_all.label(any(timelock_outcome.avg(:,timesind) > (5e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2)); % limit to right hemisphere


cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1 1]*2e-13;
cfg.channel = sensors;

figure(3); clf

cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_outcome);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('All outcomes')


cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_positive);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Positive outcome')


cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_negative);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative outcome')



cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_diff1);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
title('Negative outcome - Positive outcome')

subplot(4,2,2)
cfg = [];
cfg.xlim = times ;
cfg.zlim = [-1 1]*1e-13;
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotTFR(cfg,timelock_outcome);


subplot(4,2,4)
cfg = [];
cfg.xlim = times ;
cfg.zlim = [-1 1]*1e-13;
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotTFR(cfg,timelock_positive);


subplot(4,2,6)
cfg = [];
cfg.xlim = times;
cfg.zlim = [-1 1]*1e-13;
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotTFR(cfg,timelock_negative);

subplot(4,2,8)
cfg = [];
cfg.xlim = times ;
cfg.zlim = [-1 1]*1e-13;
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotTFR(cfg,timelock_diff1);

set(gcf,'name',['Subject ',subn(sn,:),', ',num2str(nnz(ntrials_pRPE{sn})),'/',...
    num2str(nnz(ntrials_outcome{sn})),' trials'])


%% All choices (after decision to gamble or not)
cfg = [];
cfg.trials = gamble_choices;
timelock_choice = ft_timelockanalysis(cfg,data_choice_all);
fname = sprintf('Chosen gamble option at t=0 (%d trials)',length(cfg.trials));

times = [-0.1 0.02 ; 0.02 0.13; 0.13 0.2; 0.25 0.35];
times = [-0.1 0.02 ; 0.02 0.13; 0.13 0.2; 0.25 0.4];


% find sensors with highest ER
figure(4); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_choice.time>= timee(1) & timelock_choice.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choice.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-1 1.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRT*';
timee = times(2,:);
timesind = timelock_choice.time>= timee(1) & timelock_choice.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choice.avg(:,timesind) > (5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRC*';
timee = times(3,:);
timesind = timelock_choice.time>= timee(1) & timelock_choice.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choice.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_choice);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRP*';
timee = times(4,:);
timesind = timelock_choice.time>= timee(1) & timelock_choice.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choice.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_choice);
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
    ft_topoplotTFR(cfg,timelock_choice);
end
set(gcf,'color','w','name',fname);

% All choices (after decision not to gamble)
cfg = [];
cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
% timelock_choicec = ft_timelockanalysis(cfg,data_choice_all);
fname = sprintf('Chosen certain option at t=0 (%d trials)',length(cfg.trials));
% find sensors with highest ER
figure(5); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_choicec.time>= timee(1) & timelock_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choicec.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-1 1.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRT*';
timee = times(2,:);
timesind = timelock_choicec.time>= timee(1) & timelock_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choicec.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRC*';
timee = times(3,:);
timesind = timelock_choicec.time>= timee(1) & timelock_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choicec.avg(:,timesind) < (-3e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_choicec);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)

% sensors = 'MRP*';
timee = times(4,:);
timesind = timelock_choicec.time>= timee(1) & timelock_choicec.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_choicec.avg(:,timesind) < (-4e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_choicec);
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
    ft_topoplotTFR(cfg,timelock_choicec);
end

cfg = [];
cfg.trials = 1:size(data_choice_all.time,2); cfg.trials(gamble_choices) = [];
set(gcf,'color','w','name',fname);

%% All trials (option presentation)
cfg= [];
cfg.trials = gamble_trials;
fname = sprintf('Options at t=0, chosen gamble at t in red (%d trials)',length(cfg.trials));

b = 15; % Number of bins for reaction time distribution
% timelock_trials = ft_timelockanalysis(cfg,data_trial_all); 

RT = reaction_time(gamble_trials);

times = [ 0.07 0.15; 0.18 0.25; 0.3 0.4;0.45 0.8];
% find sensors with highest ER
figure(6); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_trials.time>= timee(1) & timelock_trials.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trials.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 2.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRT*';
timee = times(2,:);
timesind = timelock_trials.time>= timee(1) & timelock_trials.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trials.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRC*';
timee = times(3,:);
timesind = timelock_trials.time>= timee(1) & timelock_trials.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trials.avg(:,timesind) < (-4e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_trials);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRP*';
timee = times(4,:);
timesind = timelock_trials.time>= timee(1) & timelock_trials.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trials.avg(:,timesind) < (-6e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_trials);
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
    ft_topoplotTFR(cfg,timelock_trials);
end
set(gcf,'color','w','name',fname);

cfg= [];
cfg.trials = 1:size(data_trial_all.time,2); cfg.trials(gamble_trials) = [];
fname = sprintf('Options at t=0, chosen certain at t in red (%d trials)',length(cfg.trials));

RT = reaction_time;
RT(gamble_trials) = [];

% timelock_trialsc = ft_timelockanalysis(cfg,data_trial_all); 
% find sensors with highest ER
figure(7); clf

% sensors ='MRF*';
timee = times(1,:);
timesind = timelock_trialsc.time>= timee(1) & timelock_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trialsc.avg(:,timesind) < (-5e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2)); % limit to right hemisphere
% sensors = 'MLT*';
cfg = [];
cfg.xlim = [-0.2 2.0];
cfg.ylim = [-1 1]*1e-13;
cfg.channel = sensors;
subplot(421)
ft_singleplotER(cfg,timelock_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRT*';
timee = times(2,:);
timesind = timelock_trialsc.time>= timee(1) & timelock_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trialsc.avg(:,timesind) > (4e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(423)
ft_singleplotER(cfg,timelock_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRC*';
timee = times(3,:);
timesind = timelock_trialsc.time>= timee(1) & timelock_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trialsc.avg(:,timesind) < (-4e-14),2 ));
sensors = sensors(strncmp(sensors,'MR',2));
cfg.channel = sensors;
subplot(425)
ft_singleplotER(cfg,timelock_trialsc);
grid on
hold on
fill([timee(1) timee(2) timee(2) timee(1)],[-1 -1 1 1],[0 1 0],'edgecolor','none','facealpha',0.5)
RTplot(RT,b)

% sensors = 'MRP*';
timee = times(4,:);
timesind = timelock_trialsc.time>= timee(1) & timelock_trialsc.time<=timee(2);
sensors = data_outcome_all.label(any(timelock_trialsc.avg(:,timesind) < (-6e-14),2 ));
sensors = sensors(strncmp(sensors,'ML',2));
cfg.channel = sensors;
subplot(427)
ft_singleplotER(cfg,timelock_trialsc);
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
    ft_topoplotTFR(cfg,timelock_trialsc);
end
set(gcf,'color','w','name',fname);
