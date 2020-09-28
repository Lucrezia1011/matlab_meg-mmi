function mmi_FC_aal_fix(data_name,twind,inducedopt,roiopt)
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid

% addpath /home/liuzzil2/fieldtrip-20190812/
% ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca

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

%% Read events

[bv_match,~] = match_triggers_fc(data_name);

% cue_match = bv_match.answer;
% choice_match = bv_match.choice;
% outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
% blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;


% inds = find(~isnan(bv.outcomeAmount));
% ind1 = inds(1)-1;
% 
% hsind = find(~isnan(bv.happySlider_response));
% mood_match.mood(hsind) =  bv.happySlider_response(hsind);
% Fsample = griddedInterpolant(find(mood_match.sample),mood_match.sample(mood_match.sample~=0),'linear');
% mood_match.sample(hsind) = Fsample(hsind);
% 
% bv = bv(inds,:);
% bv.trialNumber = (1:ntrials)'-1;

%%

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end


     
        %%
        
        twinde = twind/2*[-1,1];
        [dataout,ttdel]= define_trials(mood_match.sample(mood_match.sample~=0), data, tasktime, twinde,0);
               
        nchans = 116;
        ntrials = length(dataout.trial);
        S = str2double(sub)*ones(ntrials,1);
        
       
        %%
               
        ltvind = mood_match.bv_index(mood_match.sample~=0); % indeces start at 13
        ltvind(ttdel) = [];
 
        mood = mood_match.mood(ltvind);
        trials = mood_match.bv_index(ltvind)-12;
       
        ltvmood = table(S,trials',mood','VariableNames',...
            {'subject','trial','mood'});
                

%     save(save_name,'ltvmood',['AEC_',inducedopt{freq}],['Taec_',inducedopt{freq}]);
    save(save_name,'ltvmood','-append');
%     for freq = 1:length(inducedopt)
%         save(save_name,['AEC_',inducedopt{freq}],'-append');
%         save(save_name,['Taec_',inducedopt{freq}],'-append');
%     end
    
end

