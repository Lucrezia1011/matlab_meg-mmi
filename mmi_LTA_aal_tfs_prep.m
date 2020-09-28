function mmi_LTA_aal_tfs_prep(data_name,roiopt,filter_opt)
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


%% Co-register MRI

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

mri_name = [sub,'_anat+orig.BRIK'];

if ~exist(mri_name,'file')
    unix(['gunzip ',mri_name])
end

mri = ft_read_mri(mri_name,'dataformat','afni_brik');

tagset_shape = mri.hdr.TAGSET_NUM;
tagset_coord = mri.hdr.TAGSET_FLOATS;
tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa

tagset_p = zeros(1,3);  % Ideal orientation {RL; PA; IS}
for ii =1:3
    if strcmp(mri.hdr.Orientation(ii,:),'AP') || strcmp(mri.hdr.Orientation(ii,:),'PA')
        tagset_p(ii) = 2;
    elseif strcmp(mri.hdr.Orientation(ii,:),'LR') || strcmp(mri.hdr.Orientation(ii,:),'RL')
        tagset_p(ii) = 1;
    elseif strcmp(mri.hdr.Orientation(ii,:),'SI') || strcmp(mri.hdr.Orientation(ii,:),'IS')
        tagset_p(ii) = 3;
    end
end

m = [   -1  0   0   mri.dim(1)
    0   -1  0   mri.dim(2)
    0   0   1   1
    0   0   0   1] ;


tagset_coord = tagset_coord(tagset_p,:)'; % fiducials have shuffled coordinates

mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin

mri.transform = mri.transform/m;
fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mri = ft_volumerealign(cfg,mri);

if ~exist([sub,'_coreg.nii'],'file')
    writebrik([sub,'_coreg'],mri);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment MRI
if ~exist([processing_folder,'/headmodel.mat'],'file')
    cfg = [];
    cfg.output  = 'brain';
    segmentmri = ft_volumesegment(cfg,mri);
    
    % Head model
    
    cfg = [];
    cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, segmentmri);
    
    save([processing_folder,'/headmodel.mat'],'vol')
else
    load([processing_folder,'/headmodel.mat']);
end
sens = ft_read_sens(data_name,'senstype','meg');


%% AAL atlas
gridres = 5; % resolution of beamformer grid in mm

% Load fieldtrip 10mm MNI grid
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
template_grid = sourcemodel;
atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,sourcemodel.unit);

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

clear sourcemodel

%% Sourcemodel warp MNI grid

% sourcemodel based on 5mm grid MNI brain
cfg = [];
cfg.mri = mri;
cfg.warpmni = 'yes';
cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
cfg.unit      = 'm';
cfg.nonlinear = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);

%% Find location of AAL ROIs
R = length(sourcemodelAAL.tissuelabel);
locs = zeros(R,3);
locsAAL = cell(R,1);
for ii = 1:R
    ind = find(sourcemodelAAL.tissue == ii);
    voxc = mean(sourcemodel.pos(ind,:)); % centroid
    locs(ii,:) = voxc;
    
    locsAAL{ii} = sourcemodel.pos(ind,:);
    
end

if strcmp(roiopt,'g')
    locsc = locs;
    locs = cell2mat(locsAAL);
end
%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.sourcemodel.pos = locs;
cfg.sourcemodel.unit   = 'm';
cfg.siunits         = true;
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);

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

% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[0.5 50],filt_order,'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,50,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt

%% Beamfomer
icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

% Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value

if strcmp(roiopt,'g')
    VE = cell(1,R);
    n =0;
    for r = 1:R
        clc
        fprintf('SAM running %d/%d .\n', r,R)
        
        L = grid.leadfield( n + (1:size(locsAAL{r},1)) );
        
        VEr = zeros(data.sampleinfo(2),size(locsAAL{r},1));
        
        voxc = locsc(r,:); % centroid
        GD = zeros(1,size(locsAAL{r},1));
        for ii = 1:length(L)
            
            d = sqrt(sum((grid.pos(n+ii,:)-voxc).^2,2)); % distance from centroid
            GD(ii) = exp(-(d.^2)/1e-4); % gaussian weigthing
            lf = L{ii}; % Unit 1Am
            if GD(ii) > 0.05 && ~isempty(lf) % include voxels with weighting > 5%
                % %  G O'Neill method, equivalent to ft
                [v,d] = svd(lf'/Cr*lf);
                d = diag(d);
                jj = 2;
                
                lfo = lf*v(:,jj); % Lead field with selected orientation
                
                w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
                
                VEr(:,ii)  = GD(ii)*w'*data.trial{1};
                
            end
        end
        
        sf = corr(VEr); % check sign
        [~,ind] = max(GD);
        sf= sign(sf(ind,:));
        sf(isnan(sf)) = 0;
        VEr = VEr.*sf;
        VE{r} = sum(VEr,2);
        n = n + size(locsAAL{r},1);
    end
    
else
    L = grid.leadfield;
    VE = cell(1,R);
    for ii = 1:length(L)
        lf = L{ii}; % Unit 1Am
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        
        VE{ii}  = w'*data.trial{1};
        
        clc
        fprintf('SAM running %d/%d .\n', ii, R)
        
    end
end
%     fprintf('Beamformer finished\n' )

VE = cell2mat(VE);
% zscore VE, check dimensions
% VE = zscore(VE);
fprintf('Done.\n')

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
