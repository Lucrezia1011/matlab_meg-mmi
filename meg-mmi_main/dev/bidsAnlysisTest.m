snsaddpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

meginfo = readtable('~/MEG_participantsinfo.csv');
%%
clc
iiS = 4; % re-check 1 and 5 to 16

sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])


highpass = 0.5; 
lowpass = 300;
icaopt = 1;
plotopt = 0;
for iiN = 1:3
    data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
    if exist(data_name,'dir')
        [data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
    end
end


%% Check old processing method
data_name = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/',sdan,'MMI_mmi3_proc.ds'];

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';

data_proc = ft_preprocessing(cfg);
%%
iiN = 1;
data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
[bv_match,bv] = matchTriggers(data_name, BadSamples); % delete bad samples from recording? 


feedbackSamples = bv_match.outcome.sample(bv_match.outcome.win~=0);
bv_match.outcome.RPE(bv_match.outcome.win~=0)

figure; 
tasktime = data.time{1};
[datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);