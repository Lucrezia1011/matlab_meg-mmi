function mmiAALprep(data_name,twind)
% mmiAALprep(data_name,twind)
% Based on mmi_LTA_aal_prep
% Calculate evoked responses (lowpass 30Hz) to gamble feedback in source 
% space with AAL atlas (gaussina kernel)

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


%% Co-register MRI
gridres = 5;
mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

grid =mniLeadfields(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid

%% AAL atlas

% Load fieldtrip 10mm MNI grid
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,sourcemodel.unit);

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

clear sourcemodel

%% Find location of AAL ROIs
R = length(sourcemodelAAL.tissuelabel);
locsc = zeros(R,3); % centroid locations
locsAAL = cell(R,1); % AAL voxels
grid.inside(:) = 0;
for ii = 1:R
    ind = find(sourcemodelAAL.tissue == ii);
    grid.inside(ind) = 1;
    voxc = mean(grid.pos(ind,:)); % centroid    
    locsc(ii,:) = voxc;   
    locsAAL{ii} = grid.pos(ind,:);   
end

%%
filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

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


VE = cell(R,1);
for r = 1:R
    clc
    fprintf('SAM running %d/%d .\n', r,R)

    L = grid.leadfield(sourcemodelAAL.tissue == r);
    VEr = zeros(data.sampleinfo(2),size(locsAAL{r},1));
    
    voxc = locsc(r,:); % centroid
    GD = zeros(1,size(locsAAL{r},1));
    for ii = 1:length(L)
        
        d = sqrt(sum((locsAAL{r}(ii,:)-voxc).^2,2)); % distance from centroid
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
end

%     fprintf('Beamformer finished\n' )

VE = cell2mat(VE');
% zscore VE, check dimensions
VE = zscore(VE);
fprintf('Done.\n')

%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples); 

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
mood_match.mood(hsind-ind1) =  bv.happySlider_response(hsind);
Fsample = griddedInterpolant(find(mood_match.sample),mood_match.sample(mood_match.sample~=0),'linear');
mood_match.sample(hsind-ind1) = Fsample(hsind-ind1);

bv = bv(inds,:);
bv.trialNumber = (1:ntrials)'-1;



%% Mood
if isempty( blockmood_match)
    blockmood_match.sample = 0;
    blockmood_match.mood = 0;
end
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
Elta = Es;
Elta(2:ntrials) = EltaH(1:ntrials-1);
Rlta = A - Elta;


bestfit_name = '/data/MBDU/MEG_MMI3/data/behavioral/closed_LTA_coefs.csv';
opts = detectImportOptions(bestfit_name);
bf_pars = readtable(bestfit_name,opts); 
bestfit_sub = bf_pars(bf_pars.Var1 == str2double(sub),:);

g = bestfit_sub.gamma;
% g = .8;

sumE = zeros(ntrials,1);
sumR = zeros(ntrials,1);
sumEH = zeros(ntrials,1);
sumRH = zeros(ntrials,1);
for t = 1:ntrials
    sumEH(t) = sum( g.^(0:(t-1))' .* EltaH(t:-1:1) );
    sumRH(t) = sum( g.^(0:(t-1))' .* RltaH(t:-1:1) );
end
sumE(2:ntrials) = sumEH(1:ntrials-1);
sumR(2:ntrials) = sumRH(1:ntrials-1);

% figure; 
% plot(hsind-ind1,bv.happySlider_response(hsind-ind1),'o-')
% hold on
% plot(1:ntrials, bestfit_sub.w_LTE*sumEH +  bestfit_sub.w_LTR*sumRH + bestfit_sub.m_0)
% plot(1:ntrials, bestfit_sub.w_LTE*sumE +  bestfit_sub.w_LTR*sumR + bestfit_sub.m_0)

%%

datave = data;
datave.trial{1} = VE';
datave.label = sourcemodelAAL.tissuelabel';

% ltvall = [Elta, Rlta, sumE, sumR, Es, Rs];

fevoked = 300;

indGamble = outcome_match.win~=0 & outcome_match.sample ~=0; 
[dataout,ttdel]= define_trials(outcome_match.sample(indGamble), datave, tasktime, twind,1);

ltvind = outcome_match.bv_index(indGamble) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 300Hz for ease of memory
dataout = ft_resampledata(cfg, dataout);

nchans = length(dataout.label);
ntrials = length(dataout.trial);
S = str2double(sub)*ones(ntrials,1);

% Use Hanna's definitions
ltvcut = [EltaH(ltvind), RltaH(ltvind), sumEH(ltvind), sumRH(ltvind), Es(ltvind), Rs(ltvind)];

trials = bv.trialNumber;
trials = trials(ltvind);

xi = outcome_match.sample(indGamble);
moodi = Fmood(xi); % interpolated mood timecourse
moodi(ttdel) = [];

ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Yout = cell2mat(dataout.trial);
Yout = reshape(Yout,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);


processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
save_name = [processing_folder,'evoked_outcome_AAL'];

save(save_name ,'Yout','ltvout')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indCue = cue_match.sample~=0;
[datacue,ttdel]= define_trials(cue_match.sample(indCue), datave, tasktime, twind,1);
% Anticipatory window lasts 4s
ltvind = cue_match.bv_index(indCue) - ind1; % indeces start at 13
ltvind(ttdel) = [];
cfg = [];
cfg.resamplefs = fevoked; % Downsample to 200Hz for ease of memory
datacue = ft_resampledata(cfg, datacue);

xi = cue_match.sample(indCue);
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
ltvcut = [EltaH(ltvind-1), RltaH(ltvind-1), sumEH(ltvind-1), sumRH(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];


ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Ycue = cell2mat(datacue.trial);
Ycue = reshape(Ycue,[nchans, size(datacue.trial{1},2),length(datacue.trial)]);


processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
save_name = [processing_folder,'evoked_cue_AAL'];

save(save_name ,'Ycue','ltvcue')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indChoice = choice_match.sample~=0;
[datachoice,ttdel]= define_trials(choice_match.sample(indChoice), datave,tasktime, twind-0.3,1);

ltvind = choice_match.bv_index(indChoice) - ind1; % indeces start at 13
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
ltvcut = [EltaH(ltvind-1), RltaH(ltvind-1), sumEH(ltvind-1), sumRH(ltvind-1), Es(ltvind-1), Rs(ltvind-1)];

ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE'});

Ychoice = cell2mat(datachoice.trial);
Ychoice = reshape(Ychoice,[nchans, size(datachoice.trial{1},2),ntrials]);


processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
save_name = [processing_folder,'evoked_choice_AAL'];

save(save_name ,'Ychoice','ltvchoice')

end

