function mmiAALprep(data_name,twind, opt)
% Lucrezia Liuzzi, last updated 2021/03/15
% Based on old mmi_LTA_aal_prep
% 
% Calculate evoked responses (lowpass 30Hz) to gamble feedback, gamble 
% choice, and/or option presentation in source space with AAL atlas (gaussina kernel)
% 
% Saves timecourses and corresponding mood model parameters as .mat file
%
% mmiAALprep(data_name,twind)
% data_name = name of dataset (.ds)
% twind     = time window in seconds [t1, t2]
% opt       = cell array with trigger selection, e.g. {'outcome';'cue'}
%             'outcome': gamble feedback
%             'cue'    : gamble options presentation 
%             'choice' : gamble choice selection
% Warning: data path and output directory are hard-coded!

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

filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


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

LTAvars = LTA_calc(bv);

%% Co-register MRI
gridres = 5; % 5mm grid
mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

grid =mniLeadfields(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid

%% AAL atlas

% Load fieldtrip 5mm MNI grid
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
            
            w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ; % weights with depth correction
            
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


%% Select time windows

opt = lower(opt);

% Gamble feedback (aka outcome)
if any(strncmp(opt,'out',3)) % accept 'out' as option
    data.trial{1} = VE';
    data.label = sourcemodelAAL.tissuelabel';

    fevoked = 300; 

    indGamble = outcome_match.win~=0 & outcome_match.sample ~=0;
    [dataout,ttdel]= define_trials(outcome_match.sample(indGamble), data, tasktime, twind,1);

    ltvind = outcome_match.bv_index(indGamble) - ind1; % indeces start at 13
    ltvind(ttdel) = [];
    cfg = [];
    cfg.resamplefs = fevoked; % Downsample to 300Hz for ease of memory
    dataout = ft_resampledata(cfg, dataout);

    nchans = length(dataout.label);
    ntrials = length(dataout.trial);
    S = str2double(sub)*ones(ntrials,1);

    % Old version with E(t) sum(A_1-->t)
    ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind), LTAvars.E_sum(ltvind), ...
        LTAvars.R_sum(ltvind), LTAvars.E(ltvind), LTAvars.R(ltvind), LTAvars.M(ltvind)];

    trials = bv.trialNumber;
    trials = trials(ltvind);

    xi = outcome_match.sample(indGamble);
    moodi = Fmood(xi); % interpolated mood timecourse
    moodi(ttdel) = [];
    % save as table
    ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});

    Yout = cell2mat(dataout.trial);
    Yout = reshape(Yout,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

    save_name = [processing_folder,'evoked_outcome_',roiopt];

    save(save_name ,'Yout','ltvout')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamble option presentation
if any(strcmp(opt,'cue'))
    indCue = cue_match.sample~=0;
    [datacue,ttdel]= define_trials(cue_match.sample(indCue), data, tasktime, twind,1);
    
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
    ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind-1), LTAvars.E_sum(ltvind), ...
        LTAvars.R_sum(ltvind-1), LTAvars.E(ltvind), LTAvars.R(ltvind-1), LTAvars.M(ltvind-1)];

    ltvcue = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});

    Ycue = cell2mat(datacue.trial);
    Ycue = reshape(Ycue,[nchans, size(datacue.trial{1},2),length(datacue.trial)]);


    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
    save_name = [processing_folder,'evoked_cue_',roiopt];

    save(save_name ,'Ycue','ltvcue')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save gamble choice
if any(strcmp(opt,'choice'))
    indChoice = choice_match.sample~=0;
    % shifts time window by -0.3s
    [datachoice,ttdel]= define_trials(choice_match.sample(indChoice), data,tasktime, twind-0.3,1);

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
    ltvcut = [LTAvars.E_LTA(ltvind), LTAvars.R_LTA(ltvind-1), LTAvars.E_sum(ltvind), ...
        LTAvars.R_sum(ltvind-1), LTAvars.E(ltvind), LTAvars.R(ltvind-1), LTAvars.M(ltvind-1)];

    ltvchoice = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
        ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
        {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});

    Ychoice = cell2mat(datachoice.trial);
    Ychoice = reshape(Ychoice,[nchans, size(datachoice.trial{1},2),ntrials]);

    save_name = [processing_folder,'evoked_choice_',roiopt];

    save(save_name ,'Ychoice','ltvchoice')
end

end

