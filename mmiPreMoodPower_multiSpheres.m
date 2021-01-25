function mmiPreMoodPower_multiSpheres(data_name,roiopt,gridres,freqband,mu)
% Created August 3 2020: mmi_grid_prep_Power with new preprocessing
% roiopt = 'grid' mni grid
% gridres = grid resolution in mm
% Try with 4X weights normalization next


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

sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)

processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

highpass = 0.5;
lowpass = 300;
icaopt = 1;
plotopt = 0;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt);
f = data.fsample;
filt_order = []; % default

data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,freqband,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt

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

if strcmp(roiopt,'sens')
    
    
    [datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
    ntrials = length(datave.trial);
    %% Sensor level
    datavem = cell2mat(datave.trial);
    datas = reshape(datavem,[size(datavem,1),datave.sampleinfo(1,2),ntrials]);
    
    V = squeeze(var(datas,0,2));
    mood(ttdel) = [];
    trials(ttdel) = [];
    S = repmat(sub,length(mood),1);
    
    for iiF  = 3:7 % E,R and M from LTA model
        LTAvars.(LTAfields{iiF})(ttdel)  = [];
    end
    
    ltvmood = table(S,trials',mood',LTAvars.E_LTA ,LTAvars.E_sum,LTAvars.R_LTA ,...
        LTAvars.R_sum,LTAvars.M,'VariableNames',...
        {'subject','trial','mood','E','E_sum','RPE','RPE_sum','M'});
    
    
    save_name = sprintf('%s/pre_mood_%s_%.0f-%.0fHz',...
        processing_folder,roiopt,freqband(1),freqband(2));
    
    save(save_name,'ltvmood','V');
    
    
else
    
    %% Co-register MRI
    
    mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
    if ~exist(mri_name,'file')
        mri_name = [mri_name,'.gz'];
    end
    fids_name =  ['sub-',sub,'_fiducials.tag'];
    mri = fids2ctf(mri_name,fids_name,0);
    
    grid =mniLeadfields_multiSpheres(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid
    
    
    %%
    
    icacomps = length(data.cfg.component);
    
    C = cov(data.trial{1}');
    E = svd(C);
    nchans = length(data.label);
    noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
    % Cr = C + 4*noiseC; % old normalization, worth trying
    Cr = C + mu*E(1)*eye(size(C)); % 5% max singular value =~ 70*noise, standard
    
    [datave,ttdel]= define_trials(mood_sample, data, tasktime, [0,3],0);
    ntrials = length(datave.trial);
    Ctt = zeros([size(C),ntrials]);
    for tt=1:ntrials
        Ctt(:,:,tt) = cov(datave.trial{tt}');
    end
    
    %% Beamfomer
    
    L = grid.leadfield(grid.inside);
    W = cell(size(L));
    P = cell(size(L));
    for ii = 1:length(L)
        lf = L{ii}; % Unit 1Am
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        % depth correct
        %     w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        w = Cr\lfo / (lfo'/Cr*lfo) ;
        W{ii} = w;
        pp  = zeros(ntrials,1);
        for tt = 1:ntrials
            pp(tt) =  w'*Ctt(:,:,tt)*w;
        end
        
        P{ii} = pp/(w'*noiseC*w);
        if mod(ii,300) == 0
            clc
            fprintf('%s\n%.0f-%.0fHz: SAM running %.1f\n',...
                data_name,freqband(1),freqband(2),ii/length(L)*100)
        end
        
    end
    
    P  = cell2mat(P)';
    
    %%
    
    save_name = sprintf('%s/pre_mood_%s_%.0f-%.0fHz_mu%.0f_multiSpheres',...
        processing_folder,roiopt,freqband(1),freqband(2),mu*100);
    
    mood(ttdel) = [];
    trials(ttdel) = [];
    S = repmat(sub,length(mood),1);
    
    for iiF  = 3:7 % E,R and M from LTA model
        LTAvars.(LTAfields{iiF})(ttdel)  = [];
    end
    
    ltvmood = table(S,trials',mood',LTAvars.E_LTA ,LTAvars.E_sum,LTAvars.R_LTA ,...
        LTAvars.R_sum,LTAvars.M,'VariableNames',...
        {'subject','trial','mood','E','E_sum','RPE','RPE_sum','M'});
    %%
    save(save_name,'ltvmood','P');
    
end

%%
% W = cell2mat(W); % plkot the weights
Pgrid = zeros(size(grid.inside));
Pgrid(grid.inside) =  mean(P,2);
sourceant.pow = Pgrid;
sourceant.dim = [32 39 34]; % dimension of template
sourceant.inside = grid.inside;
sourceant.pos = grid.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'ctf';


crang = [];
cfg = [];
cfg.method        = 'slice'; %'ortho'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
ft_sourceplot(cfg, sourceout_Int);
title(data_name)
