function mmiFeedbackBF(data_name,freqband,mu,gridres,windLength,windWait,plotopt)
% mmiButtonpressBF(data_name,freqband,mu,gridres,windLength,windWait,mu,plotopt)
% For beta desynch: freqband = [13 30]; windLength = 0.5;  windWait = 0.6
% For gamma synch: freqband = [40 90]; windLength = 0.1;  windWait = -0.2
%% Co-register MRI from fiducial positions

sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)

processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];

highpass = 0.5;
lowpass = 300;
icaopt = 1;

[data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt,0);
f = data.fsample;

%% Co-register MRI

mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

grid =mniLeadfields(data_name,processing_folder,gridres,mri); % calculate leadfields on MNI grid


%% Filter data in band

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, freqband,filt_order,'but')';
data.trial{1}  =data_filt';
icacomps = length(data.cfg.component);

C = cov(data_filt);
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
Cr = C + mu*noiseC; % old normalization mu =4
% Cr = C + 0.05*E(1)*eye(size(C)); % 5% max singular value =~ 70*noise


%% Button presses 

[bv_match,~] = matchTriggers(data_name, BadSamples);

cue_match = bv_match.answer;
choice_match = bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;

indCertain = choice_match.gamble==0 & choice_match.sample~=0;
indWin = abs(outcome_match.RPE)>5 & outcome_match.sample ~=0; 


% Avoid overlap with successive sample
if windWait>0
    for tt = 1:length(samples)-1
        if samples(tt)+(windWait+windLength)*f > samples(tt+1)
            samples(tt) = 0;
        end
    end
else    
    % Avoid overlap with successive sample
    for tt = 2:length(samples)
        if samples(tt)+windWait*f < samples(tt-1)+windLength*f
            samples(tt) = 0;
        end
    end
end
samples(samples==0) = [];


[datapress,ttdel] = define_trials(outcome_match.sample(indWin),data,bv_match.time,[0 windLength],1);

Ca = zeros([size(C),length(datapress.trial)]);
for tt= 1:length(datapress.trial)
    Ca(:,:,tt) = cov(datapress.trial{tt}');
end

samples(ttdel)  = [];
[datacont,ttdel] = define_trials(outcome_match.sample(indWin),data,bv_match.time,windWait +[0 windLength],1);


Cc = zeros([size(C),length(datacont.trial)]);
for tt= 1:length(datacont.trial)
    Cc(:,:,tt) = cov(datacont.trial{tt}');
end

Ca(:,:,ttdel) = [];
Ca = mean(Ca,3);
Cc = mean(Cc,3);


%% Beamformer

W = cell(size(grid.inside));
indL = find(grid.inside)';

Tstat(1:size(grid.inside,1)) = {0};

for ii = indL
    lf = grid.leadfield{ii};
    
    % %           G O'Neill method, equivalent to ft
    [v,d] = svd(lf'/Cr*lf);
    d = diag(d);
    
    jj = 2;  % The minumum singular value is degenerate
    lfo = lf*v(:,jj); % Lead field with selected orientation
    
    w = Cr\lfo / (lfo'/Cr*lfo) ;
    %         w = Cr\lfo / sqrt((lfo'/(Cr^2)*lfo)) ;
    %         datax = w'*datafilt.trial{1};
    
    Qa = w'*Ca*w;
    Qc = w'*Cc*w;
    Tstat{ii} =(Qa - Qc) ./ (Qa + Qc); % normalised version ;
    
    %         n = w'*noiseC*w
    %         W{ii} = w;
    
    
    if mod(ii,500) == 0
        clc
        fprintf('SAM running %.1f perc.\n', ii/length(W)*100 )
    end
    
end
clc
fprintf('SAM finished\n' )
Tstat = cell2mat(Tstat);

%% Plots
if strcmp(plotopt,'anat')
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    % mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');
    
    template_grid = sourcemodel;
    clear sourcemodel
    % sourcemodel based on 5mm grid MNI brain
    cfg = [];
    cfg.mri = mri;
    cfg.warpmni = 'yes';
    cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
    cfg.unit      = 'm';
    cfg.nonlinear = 'yes';
    sourcemodel = ft_prepare_sourcemodel(cfg);
    
    
    sourceant = [];
    sourceant.pow = Tstat;
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'ctf';
    
    crang = [];
    % crang = [-.3 -.1];
    cfg = [];
    cfg.method        = 'ortho'; %'slice'
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
    
elseif strcmp(plotopt,'mni')
    %% Plot with MNI
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');
    
    sourceant = [];
    sourceant.pow = Tstat;
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'mni';
    
    
    crang = []*1e13;
    cfg = [];
    cfg.method        = 'ortho'; %'slice'
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
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
    ft_sourceplot(cfg, sourceout_Int);
end

%% Save result

save_name = sprintf('%sbuttonpressBF_%.0f-%.0fHz_a%.1fs_c%.1fs.mat',...
    processing_folder,freqband(1),freqband(2),windLength,windWait);
save(save_name,'Tstat');