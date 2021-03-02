function mmi_LTA_sensors_prep(data_name,channels)
% [1:7,11,14:16]
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions

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

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)


bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
for ii = 1:length(bv_names)
    if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
        bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
ltv = readtable(bv_name,opts); % bahavioral data


processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

bv_match = match_triggers_fc(data_name);

outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;


%% Mood

indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;

[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

xi = outcome_match.sample(outcome_match.win~=0);
moodi = Fmood(xi); % interpolated mood timecourse


%%
dataf = data;

% Evoked responses
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35,[],'but');
dataf.trial{1}= data_filt;

[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);

ltvind = outcome_match.bv_index(outcome_match.win~=0) - 12; % indeces start at 13
ltvind(ttdel) = [];
moodi(ttdel) = [];

cfg = [];
cfg.resamplefs = 200; % Downsample to 200Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);


nchans = length(dataout.label);

% Try this with all subjects!
ntrials = length(dataout.trial);
% S = sn*ones(nchans*ntrials,1);
S = str2double(sub)*ones(ntrials,1);
% sen = repmat(dataout.label,[ntrials,1]);

ltvcut = [ltv.EC, ltv.EG, ltv.Ediff, ltv.LTA, ltv.new_p, ltv.RPE, ltv.LTA_sum, ltv.RPE_sum];
ltvcut = ltvcut(ltvind,:);
trials = ltv.Var1;
trials = trials(ltvind);

% ltvcut(:,4:11) = ltvcut;
% ltvcut(:,1) = S;
% ltvcut(:,2) = trials;
% ltvcut(:,3) = moodi;

ltvall = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
        {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});


Y = cell2mat(dataout.trial);
Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

% Do this at the end??
%     % for evoked responses!! To deal with sign uncertainty
%     if ~isempty(Yall)
%         vec = corr(mean(Y,3)');
%         Y = Y.*sign(vec(:,69)); % aligns to the left precentral lobule
%         vec = corr(mean(Y,3)',mean(Yall,3)');
%         Y = Y.*sign(vec(69,69));
%     else
%         vec = corr(mean(Y,3)');
%         Y = Y.*sign(vec(:,69));
%     end
%%%%%%%%%%%%%%

% Induced oscillations
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [4 8],[],'but');
data_filt = abs(hilbert(data_filt'));
dataf.trial{1}= zscore(data_filt)';
[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
cfg = [];
cfg.resamplefs = 100; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);
Ytheta =  cell2mat(dataout.trial);
Ytheta = reshape(Ytheta,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [13 30],[],'but');
data_filt = abs(hilbert(data_filt'));
dataf.trial{1}= zscore(data_filt)';
[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
cfg = [];
cfg.resamplefs = 100; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);
Ybeta =  cell2mat(dataout.trial);
Ybeta = reshape(Ybeta,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end
save(save_name ,'Ybeta','Ytheta','Y','ltvall')

