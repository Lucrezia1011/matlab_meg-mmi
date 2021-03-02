function mmi_LTA_aal_tfs_prep_ltvonly(data_name,roiopt,filter_opt)
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
% [1:7,11,14:16]
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% Co-register MRI from fiducial positions


sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)


%% Read events

[bv_match,bv] = match_triggers_fc(data_name);

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
mood_match.mood(hsind) =  bv.happySlider_response(hsind);
Fsample = griddedInterpolant(find(mood_match.sample),mood_match.sample(mood_match.sample~=0),'linear');
mood_match.sample(hsind) = Fsample(hsind);

bv = bv(inds,:);
bv.trialNumber = (1:ntrials)'-1;


%% Mood

indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;

[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

% Standard model
A = bv.outcomeAmount; % Outcome

% Only for gamble, to check for outome window, anticioation and
% consumation
Es = (bv.winAmount + bv.loseAmount )/2;
Rs = A - Es;
risk = bv.winAmount - bv.loseAmount;
% LTA model
EltaH = cumsum(A)./(bv.trialNumber+1); % Expectation, defined by Hanna
RltaH = A - EltaH; % Assume RPE of first trial is 0
Elta = zeros(size(EltaH));
Elta(2:ntrials) = EltaH(1:ntrials-1);
Rlta = A - Elta;

g = 0.8;

sumE = zeros(ntrials,1);
sumR = zeros(ntrials,1);
sumEH = zeros(ntrials,1);
sumRH = zeros(ntrials,1);
for t = 1:ntrials
    sumE(t) = sum( g.^(0:(t-1))' .* Elta(t:-1:1) );
    sumR(t) = sum( g.^(0:(t-1))' .* Rlta(t:-1:1) );
    sumEH(t) = sum( g.^(0:(t-1))' .* EltaH(t:-1:1) );
    sumRH(t) = sum( g.^(0:(t-1))' .* RltaH(t:-1:1) );
end

% figure; subplot(121)
% plot(Es); hold all
% plot(Elta)
% plot(sumE)
% plot(EltaH)
% plot(sumEH)
% legend('Rutledge','E0','\SigmaE0','E LTA','\SigmaE LTA')
%
% subplot(122)
% plot(Rs); hold all
% plot(Rlta)
% plot(sumR)
% plot(RltaH)
% plot(sumRH)
% legend('Rutledge','R0','\SigmaR0','R LTA','\SigmaR LTA')
%%

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end

datave = data;
datave.trial{1} = VE';
datave.label = sourcemodelAAL.tissuelabel';

switch filter_opt
    case 'multi' % multitapers
        tstep = 0.05;
        cfg = [];
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.method     = 'mtmconvol';
        cfg.foi        = [1.5:0.5:4,5:14,16:2:30,35:5:55,65:5:100];
        cfg.toi        = -.2:tstep:1; % previous step of 0.05s
        cfg.keeptrials = 'yes';
        cfg.pad         = 'nextpow2';
        cfg.t_ftimwin  = 5./cfg.foi;
        cfg.tapsmofrq  = 0.3*cfg.foi;
        
        
    case 'han' % for standard hanning
        tstep = 0.025;
        cfg = [];
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.method     = 'mtmconvol';
        cfg.foi        = [1:0.5:4,5:14,16:2:50];%,45:5:150];
        cfg.toi        = -.2:tstep:1; % previous step of 0.05s
        cfg.keeptrials = 'yes';
        cfg.pad         = 'nextpow2';
        cfg.taper = 'hanning';
        cfg.t_ftimwin = ones(length(cfg.foi),1).*0.4;
end
%%%%%%%%%%%%%%% Outcome
twind = [cfg.toi(1)-2, cfg.toi(end)+2];
[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), datave, tasktime, twind,1);

TFRmult = ft_freqanalysis(cfg, dataout);
% 
% ii= 2;
% base1 = mean(squeeze(mean(TFRmult.powspctrm(:,ii,:,1:8),1)),2);
% pcolor(TFRmult.time,TFRmult.freq,(squeeze(mean(TFRmult.powspctrm(:,ii,:,:),1))-base1)./base1);
% shading interp; caxis([-0.5 0.5]);%caxis([-1 1]*1e-26);
%
% TFRhan = ft_freqanalysis(cfg, dataout);
%
%
% ii = 70;
% for tt =1:67
% base1 = mean(squeeze(mean(TFRhan.powspctrm(:,ii,:,1:8),1)),2);
% subplot(311)
% pcolor(TFRhan.time,TFRhan.freq,(squeeze(mean(TFRhan.powspctrm(tt,ii,:,:),1))-base1)./base1);
% shading interp; caxis([-2 2]);%caxis([-1 1]*5e-26);
% title(sprintf('RPE = %.2f; A = %.2f',ltvout.RPE(tt),ltvout.E(tt)+ltvout.RPE(tt)))
% subplot(312)
% plot(linspace(-0.2,1,1440),erp(ii,:,tt)); ylim([-1 1]*1e-12)
% subplot(313)
% base1 = mean(squeeze(mean(TFRmult.powspctrm(:,ii,:,1:8),1)),2);
% pcolor(TFRmult.time,TFRmult.freq,(squeeze(mean(TFRmult.powspctrm(tt,ii,:,:),1))-base1)./base1);
% shading interp; caxis([-2 2]);%caxis([-1 1]*1e-26);
%
% pause(2)
% end


ltvind = outcome_match.bv_index(outcome_match.win~=0) - ind1; % indeces start at 13
ltvind(ttdel) = [];

ntrials = length(dataout.trial);
S = str2double(sub)*ones(ntrials,1);

ltvcut = [Elta(ltvind), Rlta(ltvind), sumE(ltvind), sumR(ltvind), Es(ltvind), Rs(ltvind)];

trials = bv.trialNumber;
trials = trials(ltvind);

xi = outcome_match.sample(outcome_match.win~=0);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];

ltvout_tfs = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

tfsout = [];
tfsout.time = TFRmult.time;
tfsout.freq = TFRmult.freq;
tfsout.powspctrm = TFRmult.powspctrm;

%%%%%%%%%%%%%%%%%% Cue
twind = [cfg.toi(1)-2, cfg.toi(end)+2];

[datacue,ttdel]= define_trials(cue_match.sample(cue_match.sample~=0), datave, tasktime, twind,1);
% Anticipatory window lasts 4s

TFRmult = ft_freqanalysis(cfg, datacue);

ltvind = cue_match.bv_index(cue_match.sample~=0) - ind1; % indeces start at 13
ltvind(ttdel) = [];

xi = cue_match.sample(cue_match.sample~=0);
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
    TFRmult.powspctrm(1,:,:,:) = [];
    datacue.trial(1)=[];
end

% LTA expectation is given by sum of outcomes before current trial
ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

ntrials = length(trials);
S = str2double(sub)*ones(ntrials,1);

ltvcue_tfs = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});
tfscue = [];
tfscue.time = TFRmult.time;
tfscue.freq = TFRmult.freq;
tfscue.powspctrm = TFRmult.powspctrm;

%%%%%%%%%%%%%%%%%% Choice
cfg.toi        = -.5:tstep:0.7;
twind = [cfg.toi(1)-2, cfg.toi(end)+2];

[datachoice,ttdel]= define_trials(choice_match.sample(choice_match.sample~=0), datave,tasktime, twind,1);

TFRmult = ft_freqanalysis(cfg, datachoice);

ltvind = choice_match.bv_index(choice_match.sample~=0) - ind1; % indeces start at 13
ltvind(ttdel) = [];

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
    TFRmult.powspctrm(1,:,:,:) = [];
end

nchans = length(datachoice.label);
ntrials = length(trials);
S = str2double(sub)*ones(ntrials,1);

% Include RPE of previous trial
ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

ltvchoice_tfs = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

tfschoice = [];
tfschoice.time = TFRmult.time;
tfschoice.freq = TFRmult.freq;
tfschoice.powspctrm = TFRmult.powspctrm;

%% evoked responses

% Used to be 35Hz should make this cleaner
data_filt = ft_preproc_lowpassfilter(datave.trial{1}, datave.fsample, 25,[],'but');

datave.trial{1}= data_filt;
clear data_filt
fevoked = 300;

[dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), datave, tasktime, [-.2 1],1);

ltvind = outcome_match.bv_index(outcome_match.win~=0) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);

nchans = length(dataout.label);
ntrials = length(dataout.trial);
S = str2double(sub)*ones(ntrials,1);

ltvcut = [Elta(ltvind), Rlta(ltvind), sumE(ltvind), sumR(ltvind), Es(ltvind), Rs(ltvind)];

trials = bv.trialNumber;
trials = trials(ltvind);

xi = outcome_match.sample(outcome_match.win~=0);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];
ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Yout = cell2mat(dataout.trial);
Yout = reshape(Yout,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);


[datacue,ttdel]= define_trials(cue_match.sample(cue_match.sample~=0), datave, tasktime, [-.2 1],1);
% Anticipatory window lasts 4s
ltvind = cue_match.bv_index(cue_match.sample~=0) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 200Hz for ease of memory
datacue = ft_resampledata(cfg, datacue);

xi = cue_match.sample(cue_match.sample~=0);
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
ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Ycue = cell2mat(datacue.trial);
Ycue = reshape(Ycue,[nchans, size(datacue.trial{1},2),length(datacue.trial)]);


[datachoice,ttdel]= define_trials(choice_match.sample(choice_match.sample~=0), datave,tasktime, [-0.5 0.7],1);

ltvind = choice_match.bv_index(choice_match.sample~=0) - ind1; % indeces start at 13
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

% Include RPE of previous trial
ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Ychoice = cell2mat(datachoice.trial);
Ychoice = reshape(Ychoice,[nchans, size(datachoice.trial{1},2),ntrials]);


save(save_name ,'Yout','ltvout','Ychoice','ltvchoice','Ycue','ltvcue',...
        'tfsout','ltvout_tfs','tfschoice','ltvchoice_tfs','tfscue','ltvcue_tfs');


%%
% ii = 30;
% b = cat(1,Yout.powspctrm(:,ii,:,:),Ycue.powspctrm(:,ii,:,:) ,Ychoice.powspctrm(:,ii,:,:) );
% b = squeeze(mean(mean(b,1),4));
%
% figure; subplot(1,3,1)
% pcolor(Ycue.time, Ycue.freq, squeeze(mean(Ycue.powspctrm(:,ii,:,:),1)) - b )
% shading interp; colorbar; caxis([-1 1]*3e-27)
% subplot(1,3,2)
% pcolor(Ychoice.time, Ychoice.freq, squeeze(mean(Ychoice.powspctrm(:,ii,:,:),1))-b)
% shading interp; colorbar;  caxis([-1 1]*3e-27)
% subplot(1,3,3)
% pcolor(Yout.time, Yout.freq, squeeze(mean(Yout.powspctrm(:,ii,:,:),1)) -b)
% shading interp; colorbar;  caxis([-1 1]*3e-27)
