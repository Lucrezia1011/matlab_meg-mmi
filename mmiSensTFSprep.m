function mmiSensTFSprep(data_name,filter_opt,idx)
% mmiSensTFSprep(data_name,twind), use to do sens level TFS and evoked
% Based on mmi_LTA_aal_prep
% Calculate evoked responses (lowpass 30Hz) to gamble feedback in source 
% space with AAL atlas (gaussina kernel)

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

% filt_order = []; % default
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[0.5 100],filt_order,'but');
% data.trial{1} = data_filt;
% clear data_filt

%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples); 

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;
% 
% ntrials = nnz(~isnan(bv.outcomeAmount));
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

switch filter_opt
    case 'multi' % multitapers
        tstep = 0.025;
        cfg = [];
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.method     = 'mtmconvol';
%         cfg.foi        = [1.5:0.5:4,5:14,16:2:30,35:5:55];
        cfg.foi        = [2:55];
        cfg.toi        = -.2:tstep:1; % previous step of 0.05s
        cfg.keeptrials = 'no';
        cfg.pad         = 'nextpow2';
        cfg.t_ftimwin  = 5./cfg.foi; 
        cfg.tapsmofrq  = 0.3*cfg.foi;
        tpad = cfg.t_ftimwin(1)/2+tstep;
        twind = [cfg.toi(1)-tpad, cfg.toi(end)+tpad];
        
    case 'han' % for standard hanning
        tstep = 0.025;
        cfg = [];
        cfg.output     = 'pow';
        cfg.channel    = 'all';
        cfg.method     = 'mtmconvol';
%         cfg.foi        = [1:0.5:4,5:14,16:2:50];%,45:5:150];
        cfg.foi        = [2:55];
        cfg.toi        = -.2:tstep:1; % previous step of 0.05s
        cfg.keeptrials = 'no';
        cfg.pad         = 'nextpow2';
        cfg.taper = 'hanning';
        cfg.t_ftimwin = ones(length(cfg.foi),1).*0.4;
        tpad = cfg.t_ftimwin(1)/2+tstep;
        twind = [cfg.toi(1)-tpad, cfg.toi(end)+tpad];
end
%%%%%%%%%%%%%%% Outcome

indGamble = strcmp(bv.outcome,'win') | strcmp(bv.outcome,'lose');

if isempty(idx) % Define RPE group by subject
    R = bv.RPE(indGamble);
    A = bv.outcomeAmount(indGamble);

    indGamble = indGamble(~isnan(bv.RPE));

    % separate groups based on RPE, abs(RPE) and outcome amount
    idx = kmeans([R,A,abs(R)],3,'Replicates',30);

    s = zeros(1,3);
    for ii = 1:3
        s(ii) = sum([mean(R(idx==ii)),mean(A(idx==ii))]);
    end
    [~,idxs] = sort(s); % in increasing order
    idx(idx==idxs(1)) = -1;
    idx(idx==idxs(2)) = 0;
    idx(idx==idxs(3)) = 1;
    % figure; scatter(R(idx==1),A(idx==1))
    % hold all
    % scatter(R(idx==2),A(idx==2))
    % scatter(R(idx==3),A(idx==3))

    outcome_match.group = zeros(size(outcome_match.RPE));
    outcome_match.group(indGamble) = idx;
    outcome_match.group = outcome_match.group(1:length(outcome_match.sample));
else % Use RPE groups derived at the group level
    indGamble = indGamble(~isnan(bv.RPE));
    outcome_match.group = zeros(size(outcome_match.RPE));
    outcome_match.group(indGamble) = idx;
    outcome_match.group = outcome_match.group(1:length(outcome_match.sample));
end
indWin = outcome_match.group==1 & outcome_match.sample ~=0; 
if nnz(indWin) > 0
    [dataWin,ttdel]= define_trials(outcome_match.sample(indWin), data, tasktime, twind,1);
    TFSwin = ft_freqanalysis(cfg, dataWin);
%     b = mean(TFSwin.powspctrm,3);
%     pcolor(TFSwin.time,TFSwin.freq,(squeeze(TFSwin.powspctrm)-b')./b')
%     shading interp; caxis([-1 1]*5e-27)
    TFSwin.ntrials = length(dataWin.trial);
else
    TFSwin = [];
end


indLose = outcome_match.sample ~=0  & outcome_match.group==-1;
[dataLose,ttdel]= define_trials(outcome_match.sample(indLose), data, tasktime, twind,1);
if ~isempty(dataLose.trial)
    TFSlose = ft_freqanalysis(cfg, dataLose);
    TFSlose.ntrials = length(dataLose.trial);
else
    TFSlose = [];
end

indNeut = outcome_match.sample ~=0  & outcome_match.group==0;
if nnz(indNeut) > 0
    [dataNeut,ttdel]= define_trials(outcome_match.sample(indNeut), data, tasktime, twind,1);
    TFSneut = ft_freqanalysis(cfg, dataNeut);
    TFSneut.ntrials = length(dataNeut.trial);
else
    TFSneut = [];
end

indCertain = choice_match.gamble==0 & choice_match.sample~=0;
if nnz(indCertain) > 0
    [datacertain,ttdel] = define_trials(choice_match.sample(indCertain), data, tasktime, twind+3,1);
    cfg.toi = cfg.toi+3;
    TFScertain = ft_freqanalysis(cfg, datacertain);
    TFScertain.ntrials = length(datacertain.trial);
else
    TFScertain = [];
end

save(sprintf('%s/TFS_sens_filter-%s',processing_folder,filter_opt),...
    'TFSwin','TFSlose','TFScertain','TFSneut');

%% Plot
% TFSplot= TFSneut;
% % TFSplot.time = cat(2,TFSplot.time,TFScertain.time);
% % TFSplot.powspctrm = cat(3,TFSplot.powspctrm,TFScertain.powspctrm);
% % Equivalent to absolute baseline
% TFSplot.powspctrm = TFSplot.powspctrm - mean(TFScertain.powspctrm,3);
% 
% cfg = [];
% cfg.baseline     = 'no';%[TFScertain.time(1) TFScertain.time(end)];
% cfg.baselinetype = 'absolute';
% % cfg.channel      = data.label{ii};
% cfg.zlim         = [-1 1]*1e-27;
% cfg.xlim         = [-0.2 1];
% cfg.showlabels   = 'yes';
% cfg.layout       = 'CTF275_helmet.mat';
% cfg.interactive  = 'no';
% figure;
% ft_multiplotTFR(cfg, TFSplot); 


