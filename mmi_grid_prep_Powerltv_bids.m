function mmi_grid_prep_Powerltv_bids(data_name,freqband)
% Modification of mmi_grid_prep_Powerltv.m

% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
% roiopt = 's' sensors
% roiopt = 'grid' mni grid
% gridres = grid resolution in mm, for 'g' and 'grid' options

% addpath /home/liuzzil2/fieldtrip-20190812/
% ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca

%% Read data

sub = data_name(5:9);

cd(['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'])
highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
f = data.fsample;
    

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



%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples); % delete bad samples from recording? 

% cue_match = bv_match.answer;
% choice_match = bv_match.choice;
% outcome_match  = bv_match.outcome;
% mood_match = bv_match.ratemood;
% blockmood_match = bv_match.blockmood;
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

% Standard model
A = bv.outcomeAmount; % Outcome
A(isnan(A)) = [];
ntrials = length(A);
% LTA model
EltaH = cumsum(A)./(1:ntrials)'; % Expectation, defined by Hanna
RltaH = A - EltaH; % Assume RPE of first trial is 0

g = 0.8;

E_LTA = zeros(ntrials,1);
RPE = zeros(ntrials,1);
for t = 1:ntrials
    E_LTA(t) = sum( g.^(0:(t-1))' .* EltaH(t:-1:1) );
    RPE(t) = sum( g.^(0:(t-1))' .* RltaH(t:-1:1) );
end

E_LTA = E_LTA(trials);
RPE = RPE(trials);
EltaH = EltaH(trials);
RltaH = RltaH(trials);


%% 
[datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
ntrials = length(datave.trial);

%%

% save_name = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/%.0f-%.0fHz_%s',...
%     freqband(1),freqband(2),sub);
save_name = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/%.0f-%.0fHz_%s',...
    freqband(1),freqband(2),sub);


n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end

mood(ttdel) = [];
trials(ttdel) = [];
S = repmat(sub,length(mood),1);
RPE(ttdel) = [];
E_LTA(ttdel) = [];
EltaH(ttdel) = [];
RltaH(ttdel) = [];
ltvmood = table(S,trials',mood',EltaH,E_LTA,RltaH,RPE,'VariableNames',...
    {'subject','trial','mood','E','E_sum','RPE','RPE_sum'});

%%
save(save_name,'ltvmood','-append');

end



