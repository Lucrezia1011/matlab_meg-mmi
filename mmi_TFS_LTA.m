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
% Psub = cell(1,14);
% Rsub = cell(1,14);
npoints = 300;
Yall = [];
Ybeta = [];
Ytheta  =[];
Sall = [];
sensall = [];
RPEall = [];
ltvall = [];

for sn = [1:6,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
clearvars -except sn subn locs Yall Ybeta Ytheta Sall sensall RPEall ltvall channels
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

runs = 1; %:length(data_names)

%%
data_name = data_names{runs};

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

% if exist([data_path,'results/ft_coreg_anat.nii'],'file')
%     mri = ft_read_mri([data_path,'results/ft_coreg_anat.nii']);
%     mri.coordsys = 'ctf';
% else
% mri_name = [sub,'_anat+orig.BRIK'];
% 
% if ~exist(mri_name,'file')
%     unix(['gunzip ',mri_name])
% end
% 
% mri = ft_read_mri(mri_name,'dataformat','afni_brik');
% 
% tagset_shape = mri.hdr.TAGSET_NUM;
% tagset_coord = mri.hdr.TAGSET_FLOATS;
% tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa
% 
% tagset_p = zeros(1,3);  % Ideal orientation {RL; PA; IS}
% for ii =1:3
%     if strcmp(mri.hdr.Orientation(ii,:),'AP') || strcmp(mri.hdr.Orientation(ii,:),'PA')
%         tagset_p(ii) = 2;
%     elseif strcmp(mri.hdr.Orientation(ii,:),'LR') || strcmp(mri.hdr.Orientation(ii,:),'RL')
%         tagset_p(ii) = 1;
%     elseif strcmp(mri.hdr.Orientation(ii,:),'SI') || strcmp(mri.hdr.Orientation(ii,:),'IS')
%         tagset_p(ii) = 3;
%     end
% end
% 
% m = [   -1  0   0   mri.dim(1)
%     0   -1  0   mri.dim(2)
%     0   0   1   1
%     0   0   0   1] ;
% 
% 
% tagset_coord = tagset_coord(tagset_p,:)'; % fiducials have shuffled coordinates
% 
% mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin
% 
% mri.transform = mri.transform/m;
% fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';
% 
% cfg = [];
% cfg.method = 'fiducial';
% cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
% cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
% cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
% cfg.coordsys = 'ctf';
% cfg.viewresult = 'no';
% 
% mri = ft_volumerealign(cfg,mri);
% 
% if ~exist([sub,'_coreg.nii'],'file')
%     writebrik([sub,'_coreg'],mri);
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment MRI
% if ~exist([processing_folder,'/headmodel.mat'],'file')
%     cfg = [];
%     cfg.output  = 'brain';
%     segmentmri = ft_volumesegment(cfg,mri);
%     
%     % Head model
%     
%     cfg = [];
%     cfg.method = 'singleshell';
%     vol = ft_prepare_headmodel(cfg, segmentmri);
%     
%     save([processing_folder,'/headmodel.mat'],'vol')
% else
%     load([processing_folder,'/headmodel.mat']);
% end
% sens = ft_read_sens(data_name,'senstype','meg');


%% Calculate lead fields
% 
% 
% cfg                 = [];
% cfg.grad            = sens;
% cfg.headmodel       = vol;
% cfg.reducerank      = 2;

% cfg.channel         = {'MEG'};
% cfg.sourcemodel.pos   = locs{sn}*1e-3; 
% cfg.sourcemodel.unit   = 'm';
% cfg.siunits         = true;
% cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
% [grid] = ft_prepare_leadfield(cfg);

%% Clean data
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = channels;
% cfg.bpfreq = [110 140];
% cfg.bpfilter = 'yes';
data = ft_preprocessing(cfg);
% 
% 
% filt_order = []; % default
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [110 140],filt_order,'but');
% 
% zchans = zscore(data_filt');
% 
% for ii = 1:size(zchans,2)
%     zchans(:,ii) = smooth(abs(zchans(:,ii)),10);
% end
% 
% tempchans = false(length(data.label),1);
% for ii = 1:length(data.label)
%     tempchans(ii) = strncmp(data.label{ii},'MRT',3) || strncmp(data.label{ii},'MLT',3);
% end
% % musclez = mean(zchans(:,tempchans)>3,2);
% musclez = mean(zchans>3,2);
% 
% % figure; plot(data.time{1},musclez) %zscore >3 in temporal channels
% ind = find(musclez>0.3); %30% of all channels channels have artifacts
% muscleind = cell(1,length(ind));
% for ii = 1:length(ind)
%     muscleind{ii} = (ind(ii)-data.fsample) : (ind(ii)+data.fsample); % pad by 1s 
% end
% muscleind = unique(cell2mat(muscleind));
    
% hold on
% plot(muscleind/data.fsample,ones(1,length(muscleind))*0.4,'rx')

% 
% cfg = [];
% cfg.dataset = data_name;
% cfg.continuous = 'yes';
% % cfg.trl = cat(2, data.sampleinfo, -ones(size(data.sampleinfo,1),1)*240);
% % cfg.trl(:,1) = cfg.trl(:,1)+240;
% % cfg.trl(:,2) = cfg.trl(:,2)-240;
% cfg.trl = [241 data.sampleinfo(end)-240 -240]; 
% [cfg, artifact] = ft_artifact_muscle(cfg);
% 
% trial = 75;
% figure(1); clf 
% plot(data.sampleinfo(trial,1):data.sampleinfo(trial,2), (data.trial{trial}(98:131,:)'),'k')
% hold on
% plot(data.sampleinfo(trial,1):data.sampleinfo(trial,2), ( data.trial{trial}(227:260,:)')- 4e-13,'k')
% ind = true ;
% ii = 1;
% while ind
%     if artifact(ii,2)<=data.sampleinfo(trial,2) && artifact(ii,1)>=data.sampleinfo(trial,1) 
%     plot(artifact(ii,1):artifact(ii,2),-ones(1,(artifact(ii,2)-artifact(ii,1))+1)*2e-13, 'r','linewidth',2)    
%     elseif artifact(ii,2)>data.sampleinfo(trial,2)
%         ind = false;
%     end
%     ii=ii+1;
% end
% ylim([-8 4]*1e-13)


f = data.fsample;

if exist([processing_folder,'/ICA_artifacts.mat'],'file')
    load([processing_folder,'/ICA_artifacts.mat']);
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data_clean    = ft_rejectcomponent(cfg, comps,data);
end

data = data_clean;
clear data_clean

% 
% cfg= [];
% cfg.inwardshift = 10; % 1cm for brain headmodels
% cfg.headmodel = vol;
% cfg.template{1} = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-24071/24071MMI_mmi3_proc.ds'];    
% 
% datar = ft_megrealign(cfg, data);


%% Read events

bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
for ii = 1:length(bv_names)
    if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
        bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
ltv = readtable(bv_name,opts); % bahavioral data

bv_match = match_triggers_fc(data_name);

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;
buttonpress = bv_match.buttonpress;
tasktime = bv_match.time;
%
% pRPE = outcome_match.win == 1 ;
% nRPE = outcome_match.win == -1 ;
% pRPE_sample = outcome_match.sample(pRPE);
% nRPE_sample = outcome_match.sample(nRPE);
%%
% bvf = fieldnames(bv_match);
% 
% cellindc = [];
% for jj = 2:(length(bvf)-2) %exclude time, ITI and buttonpress
%     ind = bv_match.(bvf{jj}).sample; 
%     ind(ind==0) = [];
%     cellind = cell(1,length(ind));
%     for ii = 1:length(ind)
%         cellind{ii} = (ind(ii)-data.fsample*0.1) : (ind(ii)+data.fsample*0.5); % pad by 1s 
%     end
%     cellindc = cat(2,cellindc, cellind);
% end
% cellindc = unique(cell2mat(cellindc));
% 
% 
% ind = buttonpress;
% cellind = cell(1,length(buttonpress));
% for ii = 1:length(ind)
%     cellind{ii} = (ind(ii)-data.fsample*0.1) : (ind(ii)+data.fsample*0.5); % pad by 1s 
% end
% buttonind = unique(cell2mat(cellind));
%  
% % inddel = unique([cellindc,buttonind,muscleind]);
% % inddel(inddel<1) = [];
% % inddel(inddel>data.sampleinfo(2)) = [];
% % indkeep = 1:data.sampleinfo(end);
% % indkeep = setdiff(indkeep,inddel);
% 
% inddel = unique(muscleind); %only exclude muscle activity
% inddel(inddel<1) = [];
% inddel(inddel>data.sampleinfo(2)) = [];
% indkeep = 1:data.sampleinfo(end);
% indkeep = setdiff(indkeep,inddel);

%% Evoked fileds 
% data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 20,[],'but');
% dataf = data;
% dataf.trial{1}= data_filt; 
% 
% % Need to add which trial they represent
% [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
% ltvind = outcome_match.bv_index(outcome_match.win~=0) - 12; % indeces start at 13
% ltvind(ttdel) = [];
% 
% cfg = [];
% cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
% dataout = ft_resampledata(cfg, dataout);
% 
% % Identify Reward related positivity
% timelock_out = ft_timelockanalysis([],dataout);
% 
% [rpe,ind] = sort(outcome_match.RPE);
% n = 15;
% [datap,ttdel]= define_trials(outcome_match.sample(ind((end-n+1):end)), dataf, tasktime, [-0.5 1]);
% [datan,ttdel]= define_trials(outcome_match.sample(ind(1:n)), dataf, tasktime, [-0.5 1]);
% 
% timelock_pos = ft_timelockanalysis([],datap);
% timelock_neg = ft_timelockanalysis([],datan);
% 
% timelock_diff = timelock_pos;
% timelock_diff.avg = timelock_pos.avg - timelock_neg.avg;

%% plots
% 
% 
% figure(1); clf
% timee = [0.2 0.25];
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensors = channels(any(timelock_out.avg(:,timesind) < (-1e-13),2 ));
% cfg = [];
% cfg.xlim = [-0.5 1.0];
% cfg.ylim = [-1 1]*1.5e-13;
% cfg.channel = sensors;
% subplot(221)
% ft_singleplotER(cfg,timelock_out);
% 
% cfg = [];
% cfg.xlim = timee;
% cfg.zlim = [-1 1]*4e-14;
% cfg.layout = 'CTF275_helmet.mat';
% subplot(222)
% ft_topoplotTFR(cfg,timelock_out);
% 
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensorsd = channels(any(timelock_diff.avg(:,timesind) < (-0.2e-13),2 ));
% sensorsd = intersect(sensorsd, sensors);
% cfg = [];
% cfg.xlim = [-0.5 1.0];
% cfg.ylim = [-1 1]*1.5e-13;
% cfg.channel = sensorsd;
% subplot(223)
% ft_singleplotER(cfg,timelock_pos);
% hold on
% jj = zeros(length(sensorsd),1);
% for ii = 1:length(sensorsd)
%     jj(ii) = find(strcmp(channels,sensorsd{ii}));
% end
% plot(timelock_neg.time, mean(timelock_neg.avg(jj,:),1))
% % ft_singleplotER(cfg,timelock_neg);
% 
% 
% 
% cfg = [];
% cfg.xlim = timee;
% cfg.zlim = [-1 1]*4e-14;
% cfg.layout = 'CTF275_helmet.mat';
% subplot(224)
% ft_topoplotTFR(cfg,timelock_diff);

%%
dataf = data;
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [13 30],[],'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 20,[],'but');
% dataf.trial{1}= zscore(data_filt)'; 
dataf.trial{1}= data_filt; 


% Need to add which trial they represent
[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
ltvind = outcome_match.bv_index(outcome_match.win~=0) - 12; % indeces start at 13
ltvind(ttdel) = [];

cfg = [];
cfg.resamplefs = 200; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);

% 
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [4 8],[],'but');
% data_filt = abs(hilbert(data_filt'));
% dataf.trial{1}= zscore(data_filt)'; 
% 
% % Need to add which trial they represent
% [dataoutt,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, [-0.5 1]);
% 
% cfg = [];
% cfg.resamplefs = 200; % Downsample to 300Hz for ease of memory
% dataoutt = ft_resampledata(cfg, dataoutt);



% analysis done for all subjects together
% nchans = length(dataout.label);
% ltvcut = [ltv.EC, ltv.EG, ltv.Ediff, ltv.LTA, ltv.new_p, ltv.RPE, ltv.LTA_sum, ltv.RPE_sum];
% ltvcut = ltvcut(ltvind,:);
% 
% Y = cell2mat(dataout.trial);
% Y = reshape(Y,[nchans,length(dataout.time{1}),length(dataout.time)]);
% y = squeeze(Y(nn,tt,:));
% 
% X = table(y,ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),ltvcut(:,4),...
%         ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%         {'MEG','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
%     
% %     G = fitglm(X,'MEG ~ RPE + subject','CategoricalVars','subject');
% G = fitglm(X,'MEG ~ LTA + RPE + LTA_sum  + RPE_sum + EG + new_p');
% 
% 
% glm = fitglm(ltvcut,y)

%% Identify Reward related positivity
% timelock_out = ft_timelockanalysis([],dataout);
% 
% 
% [rpe,ind] = sort(outcome_match.RPE);
% n = 15;
% [datap,ttdel]= define_trials(outcome_match.sample(ind((end-n+1):end)), dataf, tasktime, [-0.5 1]);
% [datan,ttdel]= define_trials(outcome_match.sample(ind(1:n)), dataf, tasktime, [-0.5 1]);
% % non-reward
% [data0,ttdel]= define_trials(outcome_match.sample(ind(abs(rpe)<2 & rpe~=0)), dataf, tasktime, [-0.5 1]);
% 
% 
% timelock_pos = ft_timelockanalysis([],datap);
% timelock_neg = ft_timelockanalysis([],datan);
% timelock_neu = ft_timelockanalysis([],data0);
% 
% timelock_diff = timelock_pos;
% % timelock_diff.avg = timelock_pos.avg - timelock_neg.avg;
% timelock_diff.avg = timelock_pos.avg - timelock_neu.avg;
% 
% 
% time  = -0.4:0.2:0.6;
% figure; 
% for ii  =1:6
% subplot(3,2,ii)
% cfg = [];
% cfg.xlim = time(ii)+[0 0.2];
% cfg.ylim = [-1 1];
% cfg.zlim = [-1 1]*0.5;
% % cfg.channel = sensors;
% cfg.layout = 'CTF275_helmet.mat';
% cfg.showlabels = 'yes';
% % ft_multiplotER(cfg,timelock_out);
% ft_topoplotTFR(cfg,timelock_diff); % reward - loss
% end
% 
% figure(2); clf
% timee = [0.2 0.25];
% timesind = timelock_out.time>= timee(1) & timelock_out.time<=timee(2);
% sensors = 'MLT21';
% cfg = [];
% cfg.xlim = [-0.5 1.0];
% cfg.ylim = [-1 1]*1.5e-13;
% cfg.channel = sensors;
% subplot(321)
% ft_singleplotER(cfg,timelock_out);
% subplot(323)
% cfg.channel = 'MLT24';
% ft_singleplotER(cfg,timelock_out);
% 
% subplot(325)
% cfg.channel = 'MLP55';
% ft_singleplotER(cfg,timelock_out);
% 
% sensors = 'MRT21';
% cfg = [];
% cfg.xlim = [-0.5 1.0];
% cfg.ylim = [-1 1]*1.5e-13;
% cfg.channel = sensors;
% subplot(322)
% ft_singleplotER(cfg,timelock_out);
% subplot(324)
% cfg.channel = 'MRT24';
% ft_singleplotER(cfg,timelock_out);
% 
% subplot(326)
% cfg.channel = 'MRP55';
% ft_singleplotER(cfg,timelock_out);
%%

% figure; pcolor(dataoutm.time,1:length(dataoutm.label), dataoutm.avg)
% shading interp; colorbar
% xlim([-0.2 1]); caxis([-2 2]*1e-13)

% n=11; sen = 173;
% % figure; 
% subplot(211); plot(dataout.time{n}, dataoutm.avg(sen,:))
% title([dataout.label{sen},': Average feedback response'])
% ylim([-1 1]*3e-13)
% subplot(212); plot(dataout.time{n}, dataout.trial{n}(sen,:))
% title([dataout.label{sen},': Single feedback response'])
% ylim([-1 1]*8e-13)


nchans = length(dataout.label);

ltvcut = [ltv.EC, ltv.EG, ltv.Ediff, ltv.LTA, ltv.new_p, ltv.RPE, ltv.LTA_sum, ltv.RPE_sum];
ltvcut = ltvcut(ltvind,:);


Y =  cell2mat(dataout.trial);
Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

Ybeta = cat(3,Ybeta,Y);

% 
% Y =  cell2mat(dataoutt.trial);
% Y = reshape(Y,[nchans, size(dataoutt.trial{1},2),length(dataoutt.trial)]);
% Ytheta = cat(3,Ytheta,Y);



% Try this with all subjects!
ntrials = length(dataout.trial);
% S = sn*ones(nchans*ntrials,1);
S = sn*ones(ntrials,1);
% sen = repmat(dataout.label,[ntrials,1]);
npoints = size(dataout.trial{1},2);

Sall = cat(1,Sall,S);
ltvall = cat(1,ltvall,ltvcut);

% 
% for tt = 1:npoints
%     yt = Y(:,tt,:);
%     yt = reshape(yt,[nchans*ntrials,1]);   
%     Yall{tt} = cat(1,Yall{tt},yt);
% end
% 
% Sall = cat(1,Sall,S);
% sensall = cat(1,sensall,sen);
% RPEall = cat(1,RPEall,A);

end

%%
% Sall = num2str(Sall); % Subject should be a categorical variable

% Can downsample further when looking at oscillatory amplitude changes
% Ychan = zeros(nchans, 1.5*30,size(Ytheta,3));
% for ii = 1:size(Ytheta,3)
%     Y = Ytheta(:,:,ii)';
%     Ychan(:,:,ii) = resample(Y,30,200)';
% end

% Ychan = zeros(nchans, 1.5*60,size(Ybeta,3));
% for ii = 1:size(Ybeta,3)
%     Y = Ybeta(:,:,ii)';
%     Ychan(:,:,ii) = resample(Y,60,200)';
% end

Ychan = Ybeta;
npoints = size(Ychan,2);
glmp = cell(length(channels),npoints);
glmc = cell(length(channels),npoints);
l = size(glmp);
N = length(channels)*npoints;
% 
% changroup = {'MLC';'MLF';'MLO';'MLP';'MLT';'MRC';'MRF';'MRO';'MRP';'MRT'};
% Ychan = zeros(length(changroup),npoints,size(Yall,3));
% for n = 1:length(changroup)
%     Ychan(n,:,:) = mean(Yall(strncmp(channels,changroup{n},3),:,:),1);
% end
% 
% glm = cell(length(changroup),npoints);
% l = size(glm);
% N = l(1)*l(2);
clc
parfor n = 1:N
%     t = 0.25;
%     [~,tt] = min(abs (dataout.time{1} - t));  
    [sen,tt] = ind2sub(l,n);
%     Y = squeeze(Yall(sen,tt,:));

    Y = squeeze(Ychan(sen,tt,:));
    
%     X = table(Y,RPEall,Sall,'VariableNames',{'MEG','RPE','subject'});
    
    X = table(Y,Sall,ltvall(:,1),ltvall(:,2),ltvall(:,3),ltvall(:,4),...
        ltvall(:,5),ltvall(:,6),ltvall(:,7),ltvall(:,8),'VariableNames',...
        {'MEG','subject','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    
%     G = fitglm(X,'MEG ~ RPE + subject','CategoricalVars','subject');
    G = fitglm(X,'MEG ~ Ediff + LTA + RPE + LTA_sum  + RPE_sum + EG + new_p + subject','CategoricalVars','subject');
    
%     G = fitlme(X,'MEG ~ LTA_sum  + RPE_sum + EG + new_p +  (1|subject)');    
    
    glmp{n} =  G.Coefficients.pValue(end-6:end);
    glmc{n} =  G.Coefficients.Estimate(end-6:end);
%     glm{n} =  [G.Coefficients.Estimate(end)*1e15, G.Coefficients.pValue(end)];
    
%     lme{tt} = fitlme(tbl,'MEG~RPE+(1|subject)'); %Fixed effects for RPE 
    
%     fprintf('Done %d/%d\n',n,N)
end
fprintf('Done')

Xo = table(squeeze(Ychan(1,1,:)),Sall,ltvall(:,1),ltvall(:,2),ltvall(:,3),ltvall(:,4),...
        ltvall(:,5),ltvall(:,6),ltvall(:,7),ltvall(:,8),'VariableNames',...
        {'MEG','subject','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});


time = linspace(-0.5,1,npoints);
Gp = cell2mat(glmp);
Gp = reshape(Gp,[length(glmp{1}),size(glmp)]);

Gc = cell2mat(glmc);
Gc = reshape(Gc,[length(glmc{1}),size(glmc)]);


for x = 7:-1:1
[sen, tt ]=  find(squeeze(Gp(x,:,:))<(0.05/npoints/2));
[senu,ii,jj] = unique(sen);
figure; set(gcf,'color','w','name',Xo.Properties.VariableNames{x+3})
for sn = 1:length(senu)
%     Y = mean(Yall(senu(sn),:,:),3);
    Y = mean(Ychan(senu(sn),:,:),3);
    tx = tt(jj==sn);
    subplot(7,5,sn)
    plot(time,Y,'k')
    hold on 
    
    [g,zz] = max(abs(Gc(x,senu(sn),tx)));
   
    text(time(tx(zz))+0.05,Y(tx(zz)),num2str(Gc(x,senu(sn),tx(zz))))
    
%     tx = time(tx);
%     fill([tx(1) tx(end) tx(end) tx(1)],[-1 -1 1 1]*8e-14,[0 0 1],'facealpha',0.2,'edgecolor','none')
    plot(time(tx),Y(tx),'*r')
    title(channels(senu(sn)))
%     title(changroup(senu(sn))
end
end
%%
t = 0.25;
[~,tt] = min(abs (dataout.time{1} - t));  
Y = squeeze(Yall(:,tt,:));

% rpe = RPEall(Sall==4);
% rpe = repmat(rpe',269,1);
% y = Y(:,Sall==4);

figure; hold all
for sn = [1:6,14:16]
rpe = RPEall(Sall==sn);
y = Y(150,Sall==sn);

scatter(rpe,y)

end

% 
% 
% 
beta = zeros(1,npoints);
for tt = 1:npoints
    [~,~,betas] = fixedEffects(lme{tt});
    beta(tt) = betas.pValue(2);
%     beta(tt) = betas.tStat(2);
end

figure; semilogy(linspace(-.5,1,npoints),beta);
hold on
semilogy(linspace(-.5,1,npoints),ones(1,npoints)*0.05/npoints)

figure; semilogy(beta);


tt= 124;
yt = Yall{tt};
tbl = table(yt,sensall,RPEall,Sall,'VariableNames',{'MEG','sensor','RPE','subject'});
lmeX = fitlme(tbl,'MEG~sensor*RPE+(1|subject)');

%%


% 
% 
% nchans = length(dataout.label);
% npoints = size(dataout.trial{1},2);
% 
% RPE = ltv.RPE;
% RPE = RPE(RPE~=0);
% A = repmat(RPE(1:length(dataout.trial)),1,npoints)';
% A = A(:);
% A = repmat(A,1,nchans)';
% A = A(:);
% Y =  cell2mat(dataout.trial);
% Y = Y(:);
% 
% t = cell2mat(dataout.time);
% t = repmat(t,[nchans,1]);
% t = t(:);
% ntrials = length(dataout.trial);
% S = sn*ones(nchans*ntrials*npoints,1);
% sen = repmat(dataout.label,[ntrials*npoints,1]);
% 
% tbl = table(Y,t,sen,A,S,'VariableNames',{'MEG','time','sensor','RPE','subject'});
% %     lme{tt} = fitlme(tbl,'MEG~RPE*sensor+(sensor|subject)'); %Fixed effects for RPE and sensor
% lme = fitlme(tbl,'MEG ~ RPE*time*sensor + (sensor|subject)'); % fit for RPE, time, sensor and interaction between them, sensors vary by subject


% 
% 
% %% Beamfomer
% 
% data_cut = data;
% data_cut.time{1}(inddel) = [];
% data_cut.trial{1}(:,inddel) = [];
% data_cut.sampleinfo(2) = length(indkeep);
% 
% data_filt = ft_preproc_bandpassfilter(data_cut.trial{1}, data.fsample, [1 150],[],'but');
% data_cut.trial{1} = data_filt;
% 
% cfg = [];
% cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
% data_cut = ft_resampledata(cfg, data_cut);
% 
% % Need to adjust time to reflect cut points
% data_cut.time{1} = data.time{1};
% data_cut.time{1}(inddel) = [];
% data_cut.time{1} = downsample(data_cut.time{1},4);
% 
% icacomps = length(data.cfg.component);
% 
% C = cov(data_cut.trial{1}');
% E = svd(C);
% nchans = length(data.label);
% noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
% 
% Cr = C + 4*noiseC; % need to normalise because of ICA
% % Cr = C + 0.01*eye(nchans)*E(1);
% L = grid.leadfield;
% 
% 
% %% Power spectra for each mood rating
% [t,mood] = plot_mood(sub,false,false);
% close
% 
% %%
% nfft = 2048;
% P = cell(size(L));
% 
% tsamp = outcome_match.sample(outcome_match.sample~=0);
% tout = tasktime(tsamp);
% 
% k = 0;
% for ii = 1:length(L)
%     lf = L{ii}; % Unit 1Am
%     if ~isempty(lf)
%             k = k+1;
%             % %  G O'Neill method, equivalent to ft
%             [v,d] = svd(lf'/Cr*lf);
%             d = diag(d);
%             jj = 2;
%            
%             lfo = lf*v(:,jj); % Lead field with selected orientation
%             w = Cr\lfo / (lfo'/Cr*lfo) ;       
%             %         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
%             %         Better Hall's or normalized weights
%             wnorm = w/sqrt(w'*noiseC*w);
%             VE = wnorm'*data_cut.trial{1};   
%            
%             
% %             VE_filt = ft_preproc_bandpassfilter(VE, data_cut.fsample, [4 8],[],'but');
% %             env = abs(hilbert(VE_filt'));
% %             envd = resample(env,1,data_cut.fsample);
% %             
% %             timed = downsample(data_cut.time{1},data_cut.fsample);
% %             figure; plot(timed,envd)
%             
%             
%             Pxx = zeros(length(t{1}),nfft/2+1);
%             for tt = 1:length(t{1})
%                 tm = t{1}(tt); % mood time
%                 ind = (data_cut.time{1} > (tm-15)) & (data_cut.time{1} < (tm+15)); % 30s window around time of mood rating
% %                 tm = tasktime(outcome_match.sample)
% %                 ind = (data_cut.time{1} > (tm-15)) & (data_cut.time{1} < (tm+15));
%                 data_mood = VE(:,ind);
%                 Fxx  = fft(data_mood',nfft);
%                 Pxx(tt,:) = abs(Fxx(1:(nfft/2+1),:));
% %                 [Pxx(tt,:),Fxx] = pwelch(data_mood',[],[],nfft,data_cut.fsample);
%             end
%             P{ii} = Pxx;
%            
%     end
%     if mod(ii,300) == 0
%         clc
%         fprintf('Done %.0f perc.\n', ii/length(L)*100)
%     end
% end
% F = linspace(0,data_cut.fsample/2,nfft/2+1);
% 
% %% 
% % delete Mood_corr_*_withresp* 
% 
% k = 0;
% R = zeros(6,nnz(grid.inside));
% pvalue = zeros(6,nnz(grid.inside));
% for ii = 1:length(L)
%     lf = L{ii}; % Unit 1Am
%     if ~isempty(lf)
%         k = k+1;
%     
%         [R(1,k),pvalue(1,k)] = corr(mean(P{ii}(:,F>=1 & F<=4),2),mood{1});
%         [R(2,k),pvalue(2,k)] = corr(mean(P{ii}(:,F>=4 & F<=8),2),mood{1});
%         [R(3,k),pvalue(3,k)] = corr(mean(P{ii}(:,F>=8 & F<=13),2),mood{1});
%         [R(4,k),pvalue(4,k)] = corr(mean(P{ii}(:,F>=13 & F<=30),2),mood{1});
%         [R(5,k),pvalue(5,k)] = corr(mean(P{ii}(:,F>=30 & F<=55),2),mood{1});
%         [R(6,k),pvalue(6,k)] = corr(mean(P{ii}(:,F>=65 & F<=100),2),mood{1});
%     end
%     
% end
% 
% Psub{sn} = P;
% Rsub{sn} = R;
% 
% for sn = slist   
% %     y = cell2mat(mood_all(ii:(ii+nsets(sn)-1)));    
% %     X = comp.trial{1}(:,(a(ii)+1):a(ii+nsets(sn)));   
% %     ii = ii+nsets(sn);
% %     mdl{sn} = fitglm(X',y','linear','Distribution','gamma');      
%     y = cell2mat(mood_all(ii:(ii+nsets(sn)-1)));    
%     y = zscore(y);
%     X = comp.trial{1}(:,(a(ii)+1):a(ii+nsets(sn)));   
% %     X = VEz((a(ii)+1):a(ii+nsets(sn)),:)'; 
%     ii = ii+nsets(sn);
%     mdl{sn} = fitglm(X',y','linear','Distribution','normal'); 
% end
% 
% 
% 
% 
% % 
% % freqnames = {'delta';'theta';'alpha';'beta';'gamma_low';'gamma_high'};
% % 
% % name_nii_cat = [];
% % for k = 1:6
% %     Ranat = zeros(1,length(L));
% %     Ranat(grid.inside) = R(k,:);
% % 
% %     Panat = zeros(1,length(L));
% %     Panat(grid.inside) = pvalue(k,:);
% %     
% %     sourceant =[];
% %     sourceant.R = Ranat;
% %     sourceant.pvalue = Panat;
% %     sourceant.dim = grid.dim;
% %     sourceant.inside = grid.inside;
% %     sourceant.pos = grid.pos;
% %     cfg = [];
% %     cfg.parameter = {'R';'pvalue'};
% %     sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
% %     sourceant_Int.anatomy = sourceant_Int.R;
% % %     writebrik(['Mood_corr_',freqnames{k},'_withresp'],sourceant_Int)
% % 
% %     rname_nii = ['Mood_corr_',freqnames{k},'_withresp.nii'];
% %     ft_write_mri(rname_nii,sourceant_Int,'dataformat','nifti');
% % 
% %     sourceant_Int.anatomy = sourceant_Int.pvalue;
% % %     writebrik(['Mood_value_',freqnames{k},'_withresp'],sourceant_Int)
% %     pname_nii = ['Mood_pvalue_',freqnames{k},'_withresp.nii'];
% %     ft_write_mri(pname_nii,sourceant_Int,'dataformat','nifti');
% %     name_nii_cat = cat(2,name_nii_cat,' ',rname_nii,' ',pname_nii);
% % 
% % end
% 
% % Create sub-briks, R and pvalue alternate.
% % name_pre = 'Mood_corr_withresp.nii';
% % unix(['3dTcat -tpattern seqplus -prefix ',name_pre,' ',name_nii_cat]);
% 
% 
% end
% %%
% 
% for sn = [1,3:4,6,8:9,14]
%     sub = subn(sn,:);
%     [t,mood,trials] = plot_mood(sub,false,false);
%     n = find(trials{1} == 67);
%     if size(n,2) == 2
%         trials{1}(n(2)) = 68;
%     end
%     
%     moodsub{sn} = [trials{1}',mood{1}];
%     close
%     
%     figure(1); hold all
%     plot(trials{1},zscore(mean(Psub{sn}{2}(:,F>=30 & F<=55),2)));
%   
%     figure(2); hold all
%     plot(trials{1},zscore(mean(Psub{sn}{1}(:,F>=13 & F<=30),2)));
%     
%     figure(3); hold all
%     plot(trials{1},mood{1})
%     
% end
% 
% trials = moodsub{1}(:,1);
% mood  = zeros(size(trials));
% Pb =  zeros(size(trials));
% Pg =  zeros(size(trials));
% for t = 1:length(trials)
%     
%     moodt = [];
%     Pbeta = [];
%     Pgamma = [];
%     for sn = [1,3:4,6,8:9,14]
%         ind = moodsub{sn}(:,1) == trials(t);
%         moodt = cat(1,moodt, moodsub{sn}(ind,2));
%         
%         P = zscore(mean(Psub{sn}{1}(:,F>=13 & F<=30),2));
%         Pbeta = cat(1,Pbeta,P(ind));
%         P = zscore(mean(Psub{sn}{2}(:,F>=30 & F<=55),2));
%         Pgamma = cat(1,Pgamma,P(ind));
%     end
%     mood(t) = mean(moodt);    
%     Pb(t) = mean(Pbeta);   
%     Pg(t) = mean(Pgamma);   
%   
% end
% 
% Rb = [];
% Rg = [];
% for sn = [1,3:4,6,8:9,14]
%     Rb = cat(1,Rb,Rsub{sn}(4,1));
%     Rg = cat(1,Rg,Rsub{sn}(5,2));
% end
% 
% Rb = fisher(Rb);
% Rg = fisher(Rg);
% 
% figure(1)
% plot(trials, Pg, 'k','LineWidth' ,2)
% title(sprintf('Low gamma power timecourses: Z = %.2f (+-%.2f stdev)',mean(Rg),std(Rg)))
% ylabel('zscored low gamma power (30-55Hz)'); xlabel('trial')
% set(gcf,'color','w')
% figure(2)
% plot(trials, Pb, 'k','LineWidth' ,2)
% title(sprintf('Beta power timecourses: Z = %.2f (+-%.2f stdev)',mean(Rb),std(Rb)))
% ylabel('zscored beta power (13-30Hz)'); xlabel('trial')
% set(gcf,'color','w')
% 
% figure(3)
% plot(trials, mood, 'k','LineWidth' ,2)
% title('Mood timecourses, N = 7')
% ylabel('Mood rating'); xlabel('trial')
% set(gcf,'color','w')
% 
% 
% 
% %%
% 
% 
% % 
% % begsample = data_segmented.sampleinfo(:,1);
% % endsample = data_segmented.sampleinfo(:,2);
% % % time = ((begsample+endsample)/2) / data_segmented.fsample;
% % 
% % time = round((begsample+endsample)/2);
% % time = indkeep(time)/f;
% % figure; 
% % x = mean(P(:,:,1:4),3)';
% % pcolor(time,1:size(P,2),x)
% % shading interp; 
% % caxis([0 8]); title('delta')
% % 
% % figure; 
% % x = mean(P(:,:,4:8),3)';
% % pcolor(time,1:size(P,2),x)
% % shading interp; 
% % caxis([0 8]); title('theta')
% % 
% % figure; 
% % x = mean(P(:,:,8:13),3)';
% % pcolor(time,1:size(P,2),x)
% % shading interp; 
% % caxis([0 5]); title('alpha')
% % 
% % figure; 
% % x = mean(P(:,:,13:30),3)';
% % pcolor(time,1:size(P,2),x)
% % shading interp; 
% % caxis([0 0.5]); title('beta')
% % [t,mood] = plot_mood(sub,false,false);
% % close 
% % F = griddedInterpolant(t{1},mood{1},'pchip');
% % moodint = F(time);
% % 
% % X = [mean(P(:,:,1:4),3), mean(P(:,:,4:8),3), mean(P(:,:,8:13),3), mean(P(:,:,13:30),3)];
% % 
% % 
% % X =  mean(P(:,:,1:4),3);
% % 
% % pvalue = cell(1,size(X,2));
% % beta = cell(1,size(X,2)); 
% % parfor ii= 1:size(X,2)
% %     Xs = smooth(X(:,ii));
% %     mdl = fitglm(moodint',Xs,'linear','Distribution','normal'); 
% %     beta{ii} = mdl.Coefficients.Estimate(2);
% %     pvalue{ii} = mdl.Coefficients.pValue(2);
% % end
% % 
% % beta = cell2mat(beta);
% % pvalue = cell2mat(pvalue);
% % 
% % [psort,ind]=sort(pvalue); 
% % indsig = ind(psort<= (0.05./(length(pvalue):-1:1)));
% % psort = zeros(size(pvalue));
% % psort(indsig) = 1;
% % 
% % pvaluesig = NaN(size(grid.inside));
% % pvaluesig(grid.inside) = psort;
% % % ind = find(pvaluesig>(0.05/nnz(grid.inside)));
% % ind = find(pvaluesig~=1);
% % betasig = NaN(size(grid.inside));
% % betasig(grid.inside) = beta;
% % betasig(ind) = 0;
% % 
% % sourceant =[];
% % sourceant.pow = betasig;
% % sourceant.dim = grid.dim;
% % sourceant.inside = grid.inside;
% % sourceant.pos = grid.pos;
% % cfg = [];
% % cfg.parameter = 'pow';
% % sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
% % 
% % 
% % crang = [];
% % cfg = [];
% % cfg.method        = 'slice';
% % cfg.funparameter = 'pow';
% % cfg.maskparameter = 'pow';
% % cfg.funcolormap  = 'auto';
% % cfg.funcolorlim   = crang;
% % cfg.opacitylim = crang;
% % ft_sourceplot(cfg, sourceant_Int);
% % 
% % 
% % 
% % 
