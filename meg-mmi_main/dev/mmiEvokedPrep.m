function mmiEvokedPrep(data_name,twind,gridres,mu,fevoked)
% mmiAALprep(data_name,twind)
% roiopt = 'AAL' or 'sens'
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

filt_order = []; % default

data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,30,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


%% Read events

[bv_match,bv] = matchTriggers(data_name, BadSamples);

% cue_match = bv_match.answer;
% choice_match = bv_match.choice;
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
%     gridres = 5;
mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

grid =mniLeadfields(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid

%%


% fevoked = 50;

indGamble = outcome_match.win~=0 & outcome_match.sample ~=0;
[dataout,ttdel]= define_trials(outcome_match.sample(indGamble), data, tasktime, twind,1);

ltvind = outcome_match.bv_index(indGamble) - ind1; % indeces start at 13
ltvind(ttdel) = [];



%% Beamfomer
icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

% Cr = C + 4*noiseC; % old normalization
Cr = C + mu*eye(nchans)*E(1); % 5% max singular value

L = grid.leadfield(grid.inside);

clear data % to save disk space

ngroup = 1000;


VE = cell(ceil(length(L)/ngroup),1);
ngroupi = ngroup;
for ii = 0:ngroup:length(L)
    
    if ii + ngroup > length(L)
       ngroupi = length(L) - ii;
    end
    VEz = cell(ngroupi,1);
    parfor z= 1:ngroupi
        lf = L{ii+z}; % Unit 1Am

        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;

        lfo = lf*v(:,jj); % Lead field with selected orientation

        % depth correct
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        %         w = Cr\lfo / (lfo'/Cr*lfo) ;
        
        vez = zeros(dataout.sampleinfo(1,2),length(ltvind));
        for tt = 1:length(ltvind)
            vez(:,tt)  =(w'*dataout.trial{tt});
        end
        VEz{z} = vez;
        
    end
    
    A = cell2mat(VEz);
    A = reshape(A,dataout.sampleinfo(1,2),ngroupi,length(ltvind));
    
    dataoutVE = dataout;
    dataoutVE.label = cell(ngroupi,1);
    for  z= 1:ngroupi
        dataoutVE.label{z} = ['ve',num2str(z)];
    end
    dataoutVE.sampleinfo = repmat(dataout.sampleinfo(1,:),length(ltvind),1);
    
    dataoutVE.trial = cell(1,length(ltvind));
    for tt = 1:length(ltvind)
        dataoutVE.trial{tt} = A(:,:,tt)';
    end
    
    cfg = [];
    cfg.resamplefs = fevoked; % Downsample to 50Hz for ease of memory
    dataoutVEr = ft_resampledata(cfg, dataoutVE);
    VE{ii/ngroup+1} = dataoutVEr.trial;
    clc
    fprintf('done %.0f/%.0f\n',ii/ngroup+1,length(VE))
end

clear A
% recombine all voxels
% format: [voxels, time, trials]
VEall = [];
for z = 1:length(VE)
    A = cell2mat(VE{z});
    A = reshape(A,size(A,1),length(dataoutVEr.time{1}),length(VE{1}));
    VEall = cat(1,VEall,A);
end
fprintf('Beamformer finished\n' )

% Same problem with source orientation
A = mean(VEall,3);
[v,mii] = max(var(A,[],2));
Cc = corr(A');
VEall = VEall.*sign(Cc(:,mii));

% figure; imagesc(mean((VEall),3))
%%
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

ltvout = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),'VariableNames',...
    {'subject','trial','mood','E_LTA','RPE_LTA','E_sum','RPE_sum','E','RPE','M'});

save_name = sprintf('%s/evoked_outcome_mu-%.2f_res-%.0fmm_fs-%.0fHz.mat',...
    processing_folder,mu,gridres,fevoked);

save(save_name ,'VEall') % already have behvioral parameters saved


end

