function mmiLTAprep(data_name,twind)
% Redo LTA mesures for mmiSensprep and mmiAALprep
% mmiAALprep(data_name,twind)
% Based on mmi_LTA_aal_prep
% Calculate evoked responses (lowpass 30Hz) to gamble feedback in sensor space 

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

%% Standard pre-processing 
sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)

processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);


%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples); 

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;

ntrials = nnz(~isnan(bv.outcomeAmount));
inds = find(~isnan(bv.outcomeAmount));
ind1 = inds(1)-1;

hsind = find(~isnan(bv.happySlider_response));
mood_match.mood(hsind-ind1) =  bv.happySlider_response(hsind);
Fsample = griddedInterpolant(find(mood_match.sample),mood_match.sample(mood_match.sample~=0),'linear');
mood_match.sample(hsind-ind1) = Fsample(hsind-ind1);

bv = bv(inds,:);
bv.trialNumber = (1:ntrials)'-1;



%% Mood
if isempty( blockmood_match)
    blockmood_match.sample = 0;
    blockmood_match.mood = 0;
end
indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;

[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

LTAvars = LTA_calc(bv);
%%

fevoked = 300;

indGamble = outcome_match.win~=0 & outcome_match.sample ~=0; 
[dataout,ttdel]= define_trials(outcome_match.sample(indGamble), data, tasktime, twind,1);

ltvind = outcome_match.bv_index(indGamble) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);

nchans = length(dataout.label);
ntrials = length(dataout.trial);
S = str2double(sub)*ones(ntrials,1);

% Use Hanna's definitions
ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind), LTAvars.E_sum(ltvind), ...
    LTAvars.R_sum(ltvind), LTAvars.E(ltvind), LTAvars.R(ltvind), LTAvars.M(ltvind)];

trials = bv.trialNumber;
trials = trials(ltvind);

xi = outcome_match.sample(indGamble);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];

ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});


save_name = [processing_folder,'evoked_outcome_sens.mat'];
save(save_name ,'ltvout','-append')

save_name = [processing_folder,'evoked_outcome_AAL.mat'];
if exist(save_name,'file')
    save(save_name ,'ltvout','-append')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indCue = cue_match.sample~=0;
[datacue,ttdel]= define_trials(cue_match.sample(indCue), data, tasktime, twind,1);
% Anticipatory window lasts 4s
ltvind = cue_match.bv_index(indCue) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 200Hz for ease of memory
datacue = ft_resampledata(cfg, datacue);

xi = cue_match.sample(indCue);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];

% Skip the first choice trial: cannot test effect of RPE and LTA
% expectation
trials = bv.trialNumber;
trials = trials(ltvind);
if trials(1) == 0
    ltvind(1) = [];
    trials(1) = [];
    moodi(1) = [];
    datacue.trial(1)=[];
end

nchans = length(datacue.label);
ntrials = length(trials);
S = str2double(sub)*ones(ntrials,1);

% Include RPE and expectation of previous trial
ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind-1), LTAvars.E_sum(ltvind), ...
    LTAvars.R_sum(ltvind-1), LTAvars.E(ltvind), LTAvars.R(ltvind-1), LTAvars.M(ltvind-1)];


ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});


save_name = [processing_folder,'evoked_cue_sens.mat'];

save(save_name ,'ltvcue','-append')

save_name = [processing_folder,'evoked_cue_AAL.mat'];
if exist(save_name,'file')
    save(save_name ,'ltvcue','-append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indChoice = choice_match.sample~=0;
[datachoice,ttdel]= define_trials(choice_match.sample(indChoice), data,tasktime, twind-0.3,1);

ltvind = choice_match.bv_index(indChoice) - ind1; % indeces start at 13
ltvind(ttdel) = [];

cfg = [];
cfg.resamplefs = fevoked; % Downsample to 200Hz for ease of memory
datachoice = ft_resampledata(cfg, datachoice);


xi = choice_match.sample(choice_match.sample~=0);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];

% Skip the first choice trial: cannot test effect of RPE and LTA
% expectation
trials = bv.trialNumber;
trials = trials(ltvind);
if trials(1) == 0
    ltvind(1) = [];
    trials(1) = [];
    moodi(1) = [];
    datachoice.trial(1)=[];
end

nchans = length(datachoice.label);
ntrials = length(trials);
S = str2double(sub)*ones(ntrials,1);

ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind-1), LTAvars.E_sum(ltvind), ...
    LTAvars.R_sum(ltvind-1), LTAvars.E(ltvind), LTAvars.R(ltvind-1), LTAvars.M(ltvind-1)];

ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});


save_name = [processing_folder,'evoked_choice_sens.mat'];

save(save_name ,'ltvchoice','-append')

save_name = [processing_folder,'evoked_choice_AAL.mat'];
if exist(save_name,'file')
    save(save_name ,'ltvchoice','-append')
end

%%  PreMoodPower
[bv_match,bv] = matchTriggers(data_name, BadSamples); 

tasktime = bv_match.time;

mood_sample = bv_match.ratemood.sample(bv_match.ratemood.sample~=0);
% mood_sample = cat(2,mood_sample,bv_match.blockmood.sample(bv_match.blockmood.sample~=0));

[mood_sample, moodind] = sort(mood_sample);

mood =  bv_match.ratemood.mood(bv_match.ratemood.sample~=0);
% mood = cat(2,mood,bv_match.blockmood.mood(bv_match.blockmood.sample~=0));

mood = mood(moodind);

trials =  bv_match.ratemood.bv_index(bv_match.ratemood.sample~=0);
% trials = cat(2,trials,bv_match.blockmood.bv_index(bv_match.blockmood.sample~=0)-0.5);

trials = trials(moodind)-12;

LTAvars = LTA_calc(bv);
LTAfields = fieldnames(LTAvars,'-full');

for iiF  = 3:7 % E,R and M from LTA model
    LTAvars.(LTAfields{iiF})  = LTAvars.(LTAfields{iiF})(trials);
end

[datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
ntrials = length(datave.trial);

mood(ttdel) = [];
trials(ttdel) = [];
S = repmat(sub,length(mood),1);

for iiF  = 3:7 % E,R and M from LTA model
    LTAvars.(LTAfields{iiF})(ttdel)  = [];
end

ltvmood = table(S,trials',mood',LTAvars.E_LTA ,LTAvars.E_sum,LTAvars.R_LTA ,...
    LTAvars.R_sum,LTAvars.M,'VariableNames',...
    {'subject','trial','mood','E','E_sum','RPE','RPE_sum','M'});

roiopt = {'sens';'grid'};
for ii = 1:2
    save_name = sprintf('%s/pre_mood_%s_%.0f-%.0fHz.mat',...
        processing_folder,roiopt{ii},25,40);
    if exist(save_name,'file')
        save(save_name,'ltvmood','-append');
    end
end


%% P300

[bv_match,bv] = matchTriggers(data_name, BadSamples); 

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
if isempty(blockmood_match)
    blockmood_match.sample = [];
    blockmood_match.mood = [];
end
tasktime = bv_match.time;


indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;
[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

% xi = gamble_sample;
% mood = Fmood(xi); % interpolated mood timecourse
%   
% trials =  choice_match.bv_index(choice_match.gamble==1)-12;

LTAvars = LTA_calc(bv);
LTAfields = fieldnames(LTAvars,'-full');

cue_sample = cue_match.sample(cue_match.choice~=0);
[datave,ttdel]= define_trials(cue_sample, data, tasktime, [.25,.4],0);
ntrials = length(datave.trial);

trials =  cue_match.bv_index(cue_match.choice~=0)-12;

% do mood and expectation (all 3 types) influence the P300?
% does not make too much sense to include RPE
for iiF  = 1:7 % E,R and M from LTA model
    LTAvars.(LTAfields{iiF}) = LTAvars.(LTAfields{iiF})(trials);
    LTAvars.(LTAfields{iiF})(ttdel)  =[];
end

trials_cue = trials;
trials_cue(ttdel) = [];

xi = cue_sample;
mood_cue = Fmood(xi); % interpolated mood timecourse
mood_cue(ttdel) = [];
Scue = repmat(sub,length(mood_cue),1);
gamble_cue = cue_match.choice(cue_match.choice~=0);
gamble_cue(ttdel) = [];

outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
for ii = 1:2
save_name = sprintf('%s/cueP300_%s.mat',outpath,roiopt{ii});

if exist(save_name,'file')

ltvcue = table(Scue,trials_cue',mood_cue',LTAvars.E ,LTAvars.E_LTA ,...
    LTAvars.E_sum,gamble_cue',LTAvars.M,...
    'VariableNames',{'subject','trial','mood','E','E_LTA','E_sum','choice','M'});
save(save_name,'ltvcue','-append');
end
end
