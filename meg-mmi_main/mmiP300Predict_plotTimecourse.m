% Made August 4th 2020
% Based on mmi_P300_predict.m

clear all
close all
clc

data_list = [];

meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
data_list = [];

% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
subexclude = [10,24];

analy_case = 'confirm';
roiopt = 'sens'; % running for grid
switch roiopt
    case 'grid'
        subexclude = [subexclude,26,49,53];
end

switch analy_case
    case 'confirm'
        subexclude = [subexclude,1:12,14:16];
    case 'explore'
        subexclude = [subexclude,13,17:56];
end

Nlist(subexclude) = []; 
zz= 0;

for sn = Nlist % all subjects with continuos recordings and latent variables
        
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

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% For sensor based analysis
load /data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat
sensall = channels;

load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/confirm/sig_channels_E_LTA.mat');
sensp = C;
load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/confirm/sig_channels_E_sum.mat');
sensp = unique([sensp;C]);
%% For grid based analysis
gridres= 5;

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/';
Y = dlmread([data_path,'/meg_trials_evoked_cue.txt']);
Y = reshape(Y, 360,116,size(Y,2));

opts = detectImportOptions([data_path,'latent_vars_evoked_outcome.csv']);
Xv = readtable([data_path,'latent_vars_evoked_outcome.csv'],opts);

%% It's impossible to do better than this for source sign. The sensor level analysis is the confirmation
for s = 1:length(data_list)
    
%% Co-register MRI from fiducial positions

data_name = data_list{s};
sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)
processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
save_name = sprintf('%s/cueP300_sens_Timecourse.mat',outpath);
if ~exist(save_name,'file')

highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
f = data.fsample;


filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


%% Get source localized data on the AAL for this subject
subs = unique(Xv.subject);
Ym = Y(:,:,subs == str2double(sub));

%% Read events

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

%%
[~,~,ind]= intersect(sensp,data.label);
% include all available trials (even ones when no choice was made)
cue_sample = cue_match.sample(cue_match.choice~=0);
[datave,ttdel]= define_trials(cue_sample, data, tasktime, [-.5,2],0);
ntrials = length(datave.trial);
datavem = cell2mat(datave.trial);

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

datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
sens_cue = mean(datas(ind,:,:),3);

save(save_name,'sens_cue','ntrials')
end
end
return

%% Plot
time = linspace(-.5,2,1200*2.5);
data = zeros(20,1200*2.5,length(Nlist));
s = 0;
for sn = Nlist % all subjects with continuos recordings and latent variables
    
    sub = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sub,'/meg'])
    
    data_list = [];
    zz = 0;
    ntot =0;
    datas = [];
    for iiN = 1:3 
        data_name = ['sub-',sub,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;                    
            outpath = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
            save_name = sprintf('%s/cueP300_sens_Timecourse.mat',outpath);
            load(save_name)
            datas{zz} = sens_cue*ntrials;
            ntot = ntot + ntrials;

        end
    end
    datas = cell2mat(datas);
    datas = reshape(datas, [size(datas,1),length(time),zz]);
    datas = sum(datas,3)/ntot;
    data(:,:,find(sn == Nlist)) = datas;
end

datal = mean(data(strncmp(sensp,'ML',2),:,:));
dataml = mean(datal,3);
datasdl = std(datal,0,3);

datar = mean(data(strncmp(sensp,'MR',2),:,:));
datamr = mean(datar,3);
datasdr = std(datar,0,3);

%%
figure; set(gcf,'color','w')
plot(time,dataml*1e15)
hold on
plot(time,-datamr*1e15)
fill([0.25 0.4 .4 0.25],[-1 -1 1 1]*1e2,[0.5 0.5 0.5],'facealpha',0.3,'edgecolor','none')
fill([time fliplr(time)],[dataml-datasdl/sqrt(nnz(Nlist)), ...
    fliplr(dataml+datasdl/sqrt(nnz(Nlist)))]*1e15,[0 0 1],'facealpha',0.3,'edgecolor','none')
fill([time fliplr(time)],-[datamr-datasdr/sqrt(nnz(Nlist)), ...
    fliplr(datamr+datasdr/sqrt(nnz(Nlist)))]*1e15,[1 0 0],'facealpha',0.3,'edgecolor','none')

grid on
ylim([-1 1]*0.6e2); xlim([-.2 ,1])
xlabel('time (s)'); ylabel('Magnetic Field (fT)')
title(sprintf('Average signal (no. subjects=%d) from significant sensors',nnz(Nlist)))
legend(['ML sensors n=',num2str(nnz(strncmp(sensp,'ML',2)))],...
    ['MR sensors (sign flipped) n=',num2str(nnz(strncmp(sensp,'MR',2)))],...
    'Hypothesis time window',...
    'ML standard error over subjects','MR standard error over subjects',...
    'location','best') % standard error over subjects

saveas(gcf,sprintf('~/matlab/figures/P300_sens_timecourse.png'))
