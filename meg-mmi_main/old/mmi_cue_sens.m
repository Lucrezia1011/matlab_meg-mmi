
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

data_name = '24071MMI_mmi3_proc.ds';
data_name = '23490MMI_mmi3_proc.ds';

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming/'];
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

if exist([processing_folder,'/ICA_artifacts.mat'],'file')
    load([processing_folder,'/ICA_artifacts.mat']);
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data          = ft_rejectcomponent(cfg, comps,data);

%%
filt_order = []; % default

% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[0.5 40],filt_order,'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


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


%%
twind = [1,4];
[dataexp,ttdel]= define_trials(choice_match.sample(choice_match.gamble==1), data, tasktime, twind,1);

ltvind = choice_match.bv_index(choice_match.gamble==1) - ind1; % indeces start at 13
ltvind(ttdel) = [];

nchans = length(dataexp.label);
ntrials = length(dataexp.trial);

for tt = 1:ntrials
    dataexp.trial{tt} = abs(dataexp.trial{tt});
end

trials = bv.trialNumber;
trials = trials(ltvind);

trials1 = find(ltvind <= 27);
trials2 = find(ltvind>27 & ltvind <= 54);
trials3 = find(ltvind>54);

nt = min([nnz(trials1),nnz(trials2),nnz(trials3)]);
trials1 = trials1(end-nt+1:end);
trials2 = trials2(end-nt+1:end);
trials3 = trials3(end-nt+1:end);

% figure; plot(1:81, sumE)
% hold on
% plot(1:81,Elta)

figure(2); clf
for tt = 1:3
cfg = [];
switch tt
    case 1
        cfg.trials = trials1;
    case 2
        cfg.trials = trials2;
    case 3
        cfg.trials = trials3;
end

timelock = ft_timelockanalysis(cfg,dataexp);

% figure; pcolor(linspace(1,4,f*3),1:271,timelock.avg)
% shading flat
subplot(1,3,tt)
cfg = [];
cfg.xlim = [];
cfg.zlim = [0.3 1.4]*1e-13;
cfg.channel = 'MEG';
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotER(cfg,timelock);
end

    % Anticipatory window lasts 4s
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

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end



% Evoked responses
if istrue(evokedopt)
   
    % Used to be 35Hz should make this cleaner
    data_filt = ft_preproc_lowpassfilter(datave.trial{1}, data.fsample, 20,[],'but');
    
    dataf.trial{1}= data_filt;

    fevoked = 300;
    
    [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, twind,1);
  
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
%     ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%             {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
            {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

    Yout = cell2mat(dataout.trial);
    Yout = reshape(Yout,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

    
    [datacue,ttdel]= define_trials(cue_match.sample(cue_match.sample~=0), data, tasktime, twind,1);
   
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
    ltvcut = [Elta(ltvind-1), Rlta(ltvind-1), sumE(ltvind-1), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

    
%     ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%         {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});
    
    Ycue = cell2mat(datacue.trial);
    Ycue = reshape(Ycue,[nchans, size(datacue.trial{1},2),length(datacue.trial)]);
    
    
    [datachoice,ttdel]= define_trials(choice_match.sample(choice_match.sample~=0), dataf,tasktime, [-0.5 0.7],1);
    
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
%     ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind), Rs(ltvind-1)];
    ltvcut = [Elta(ltvind-1), Rlta(ltvind-1), sumE(ltvind-1), sumR(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

    
%     ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%         {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});    
    
    Ychoice = cell2mat(datachoice.trial);
    Ychoice = reshape(Ychoice,[nchans, size(datachoice.trial{1},2),ntrials]);
    
    save(save_name ,'Yout','ltvout','Ychoice','ltvchoice','Ycue','ltvcue')

% Check how similar the gaussian decay is to the centroid
% n =32;
% figure(2); clf
% plot(linspace(twind(1),twind(2),1200),zscore(mean(Y(n,:,:),3)))
% hold on
% plot(linspace(twind(1),twind(2),1200),zscore(mean(Y1(n,:,:),3)))
% legend('centroid','gaussian decay')
% title(sourcemodelAAL.tissuelabel{n})
% xlim([-.2 1])
% 
% plot(abs(diag(corr(mean(Y,3)',mean(Y1,3)'))))


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
end


if ~isempty(inducedopt)
   
    clear ltvchoice ltvout ltvcue Ychoice Yout Ycue Yalpha Ybeta Ytheta 
    
    for freq = 1:length(inducedopt)
    % Induced oscillations
    switch inducedopt{freq}
        case 'delta'
            freql = [1 4];
            filtertype = 'firls';
        case 'theta'
            freql = [4 8];
            filtertype = 'but';
         case 'alpha'
            freql = [8 13];
            filtertype = 'but';
        case 'beta'
            freql = [13 30];
            filtertype = 'but';
    end
   
    finduced  = 60;
    data_filt = ft_preproc_bandpassfilter(datave.trial{1}, datave.fsample, freql,[],filtertype);
    data_filt = abs(hilbert(data_filt'));
    dataf.trial{1}= zscore(data_filt)';
    
    %%
    twind = [-3,3];
%     twindf = twind*f./mean(freql);
    
    [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), dataf, tasktime, twind,1);
        
    ltvind = outcome_match.bv_index(outcome_match.win~=0) - ind1; % indeces start at 13
    ltvind(ttdel) = [];
    cfg = [];
    cfg.resamplefs = finduced; % Downsample to 200Hz for ease of memory
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
%     ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%             {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

        
    Yout = cell2mat(dataout.trial);
    Yout = reshape(Yout,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

%     twind = [-3,3];
%     [datacue,ttdel]= define_trials(cue_match.sample(cue_match.sample~=0), dataf, tasktime, twind,0);
%     % Anticipatory window lasts 4s
%     ltvind = cue_match.bv_index(cue_match.sample~=0) - ind1; % indeces start at 13
%     ltvind(ttdel) = [];
%     cfg = [];
%     cfg.resamplefs = finduced; % Downsample to 200Hz for ease of memory
%     datacue = ft_resampledata(cfg, datacue);
% 
%     nchans = length(datacue.label);
%     ntrials = length(datacue.trial);
%     S = str2double(sub)*ones(ntrials,1);
% 
%     ltvcut = ltvall(ltvind,:);
%     trials = bv.trialNumber;
%     trials = trials(ltvind);
% 
%     xi = cue_match.sample(cue_match.sample~=0);
%     moodi = Fmood(xi); % interpolated mood timecourse
%     moodi(ttdel) = [];
%     ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
%         {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});
% 
%     Ycue = cell2mat(datacue.trial);
%     Ycue = reshape(Ycue,[nchans, size(datacue.trial{1},2),length(datacue.trial)]);
    
    twind = [-3,3];
    [datachoice,ttdel]= define_trials(choice_match.sample(choice_match.sample~=0), dataf, tasktime, twind,1);
    
    ltvind = choice_match.bv_index(choice_match.sample~=0) - ind1; % indeces start at 13
    ltvind(ttdel) = [];
    cfg = [];
    cfg.resamplefs = finduced; % Downsample to 200Hz for ease of memory
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
    ltvcut = [Elta(ltvind), Rlta(ltvind-1), sumE(ltvind), sumR(ltvind-1), Es(ltvind), Rs(ltvind-1)];
    
%     ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
%         ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
%             {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});
    ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

    Ychoice = cell2mat(datachoice.trial);
    Ychoice = reshape(Ychoice,[nchans, size(datachoice.trial{1},2),length(datachoice.trial)]);
    %%
%     figure(freq); clf; set(gcf,'color','w'); 
%     subplot(131)
%     pcolor(linspace(-2.5,2.5,100*5),1:116,mean(Ycue,3))
%     shading interp; colorbar; caxis([-.5 .5])
%     xlabel('Time (s)'); ylabel('ROI'); title('Cue')
%     subplot(132)
%     pcolor(linspace(-1,3,100*4),1:116,mean(Ychoice,3))
%     shading interp; colorbar; caxis([-.5 .5])
%     xlabel('Time (s)'); ylabel('ROI'); title('Choice')
%     subplot(133)
%     pcolor(linspace(-4,3,100*7),1:116,mean(Yout,3))
%     shading interp; colorbar; caxis([-.5 .5])
%     xlabel('Time (s)'); ylabel('ROI'); title('Outcome')
    %%
    eval(sprintf('Ychoice_%s = Ychoice; ltvchoice_%s = ltvchoice;',inducedopt{freq},inducedopt{freq}));
    eval(sprintf('Yout_%s = Yout; ltvout_%s = ltvout;',inducedopt{freq},inducedopt{freq}));
%     eval(sprintf('Ycue_%s = Ycue; ltvcue_%s = ltvcue;',inducedopt{freq},inducedopt{freq}));
    
    clear ltvchoice ltvout ltvcue Ychoice Yout Ycue
    
%     save(save_name,'-append',...
%         ['Ycue_',inducedopt{freq}],['Ychoice_',inducedopt{freq}],['Yout_',inducedopt{freq}],...
%         ['ltvcue_',inducedopt{freq}],['ltvchoice_',inducedopt{freq}],['ltvout_',inducedopt{freq}]);
    save(save_name,'-append',...
        ['Ychoice_',inducedopt{freq}],['Yout_',inducedopt{freq}],...
        ['ltvchoice_',inducedopt{freq}],['ltvout_',inducedopt{freq}]);

    end

end

