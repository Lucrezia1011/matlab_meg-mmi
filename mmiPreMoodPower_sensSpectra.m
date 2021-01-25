clear all
close all
% clc
meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

data_list = [];


% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
subexclude = [10,24];

roiopt = 'sens'; % running for grid
switch roiopt
    case 'grid'
        subexclude = [subexclude,26,49,53];
end

Nlist(subexclude) = []; 
zz= 0;
for sn = Nlist %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sdan = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;
            data_list{zz} = data_name;
        end
    end
    
   
end

%%
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat');

% MRC52, MLC16
Ps = cell(length(data_list),2);

for ii = 1:length(data_list)
    data_name = data_list{ii};
    sub = data_name(5:9);

% Created January 22 2021: 
% Caluculate power spectrum from sensors during pre-mood

data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)

processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
f = data.fsample;


%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples);
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

LTAvars = LTA_calc(bv);
LTAfields = fieldnames(LTAvars,'-full');

for iiF  = 3:7 % E,R and M from LTA model
    LTAvars.(LTAfields{iiF})  = LTAvars.(LTAfields{iiF})(trials);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
ntrials = length(datave.trial);
%% Sensor level
% ssen = strcmp(datave.label,'MRC52');
% Average over all central sensors
ssen = strncmp(datave.label,'MRC',3) | strncmp(datave.label,'MLC',3) ;

P = cell(1,ntrials);
Pf = cell(1,ntrials);
for iiS = 1:ntrials
    [P{iiS},ff] = pwelch(datave.trial{iiS}(ssen,:)',[],[],[],f);
    Pf{iiS} = fft(datave.trial{iiS}(ssen,:)');
end

P = cell2mat(P);
P = reshape(P,[size(P,1),nnz(ssen),ntrials]);
P = squeeze(mean(P,2));
Pf = abs(cell2mat(Pf));
Pf = reshape(Pf,[size(Pf,1),nnz(ssen),ntrials]);
Pf = squeeze(mean(Pf,2));
Pf = Pf(1:(size(Pf,1)/2+1),:);
fff = linspace(0,f,size(Pf,1))';

Ps{ii,1}  = P;
Ps{ii,2} = Pf;
%%
% figure(1); clf;
% subplot(211)
% plot(ff,mean(P,2))
% 
% hold on
% fill([ff;flipud(ff)],[mean(P,2)+std(P,0,2); flipud(mean(P,2)-std(P,0,2))],...
%     [0 0 1],'facealpha',0.3,'Edgecolor','none')
% xlim([0 55]); grid on
% 
% subplot(212)
% plot(fff,mean(mean(Pf,2),3))
% 
% hold on
% fill([fff;flipud(fff)],[mean(Pf,2)+std(Pf,0,2); flipud(mean(Pf,2)-std(Pf,0,2))],...
%     [0 0 1],'facealpha',0.3,'Edgecolor','none')
% xlim([0 55]); grid on
%%


end

%%
PP = zeros(size(Ps{1,1},1),length(Nlist));
PPf = zeros(size(Ps{1,2},1),length(Nlist));

for sn = Nlist %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
    
    sdan = num2str(meginfo.SDAN(sn));
    n = strncmp(data_list,['sub-',sdan],9);
    p1 = [];
    p2 = [];
    for nn = find(n)
        p1 = cat(2,p1,Ps{nn,1});
        p2 = cat(2,p2,Ps{nn,2});
    end
    PP(:,sn) = mean(p1,2);
    PPf(:,sn) = mean(p2,2);
    
end


figure(1); clf;
subplot(211)
plot(ff,mean(PP,2))

hold on
fill([ff;flipud(ff)],[mean(PP,2)+std(PP,0,2); flipud(mean(PP,2)-std(PP,0,2))],...
    [0 0 1],'facealpha',0.3,'Edgecolor','none')
xlim([0 55]); grid on

subplot(212)
plot(fff,mean(mean(PPf,2),3))

hold on
fill([fff;flipud(fff)],[mean(PPf,2)+std(PPf,0,2); flipud(mean(PPf,2)-std(PPf,0,2))],...
    [0 0 1],'facealpha',0.3,'Edgecolor','none')
xlim([0 55]); grid on
%