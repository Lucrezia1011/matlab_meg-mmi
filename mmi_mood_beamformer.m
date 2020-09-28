function sourcenorm = mmi_mood_beamformer(sub,freqband,mu,dataprep)


data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

%% Read events
bv_match = dataprep.bv;

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;

pRPE = outcome_match.win == 1 ;
nRPE = outcome_match.win == -1 ;
pRPE_sample = outcome_match.sample(pRPE);
nRPE_sample = outcome_match.sample(nRPE);

%% Filter data in theta band
% freqband = [4,8];
data_clean = dataprep.data;
f = data_clean.fsample;
filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1},f, freqband,filt_order,'but')';

C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));

% Cr = C + mu*noise;
Cr = C + mu*max(svd(C)).*eye(size(C));

%% Cue-P3 ERP: 200ms
%
% wind = 0.2*f:0.4*f;
% % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
%
% % win_logic = bv.RPE > 10 & bv.outcomeAmount > 10;
% win_logic = bv.RPE > 0;
% win_logic = win_logic(1:length(outcome_match.sample));
% win_samples = outcome_match.sample(win_logic);
%
% % lose_logic = bv.RPE < -10 & bv.outcomeAmount < -10;
% lose_logic = bv.RPE < 0;
% lose_logic = lose_logic(1:length(outcome_match.sample));
% lose_samples = outcome_match.sample(lose_logic);
%
% Cap = zeros([size(C),size(win_samples,2)]);
% for tt= 1:size(win_samples,2)
%     data_win = data_filt(win_samples(tt)+wind,:);
%     Cap(:,:,ii) = cov(data_win);
% end
% Cap = mean(Cap,3);=
%
% Can = zeros([size(C),size(lose_samples,2)]);
% for tt= 1:size(lose _samples,2)
%     data_lose = data_filt(lose_samples(tt)+wind,:);
%     Can(:,:,ii) = cov(data_lose);
% end
% Can = mean(Can,3);
%


%% Mood decision 3s

ind = mood_match.sample>0;
mood_match.mood(~ind)= NaN;

v1 = mood_match.mood(ind);
samples1 = mood_match.sample(ind);

ind = blockmood_match.sample>0;
blockmood_match.mood(~ind) = NaN;
v0 = blockmood_match.mood(ind);
samples0 = blockmood_match.sample(ind);

v = [v1,v0];
[s,ind] = sort(v);

samples = [samples1, samples0];
samples = samples(ind);

hs = s > median(s);
ls = s < median(s);

happy_samples = samples(hs);
sad_samples = samples(ls);

if nnz(hs)> nnz(ls)
    happy_samples(1) = [];
elseif nnz(hs)< nnz(ls)
    sad_samples(end) = [];
end

wind = 0*f:3*f; % Exclude edges for sensory and motor overlap
% cuep3_sample = zeros(nnz(choice_match.sample), length(wind));


Cah = zeros([size(C),size(happy_samples,2)]);
for tt= 1:size(happy_samples,2)
    data_win = data_filt(happy_samples(tt)+wind,:);
    Cah(:,:,tt) = cov(data_win);
end
Cah = mean(Cah,3);

Cal = zeros([size(C),size(sad_samples,2)]);
for tt= 1:size(sad_samples,2)
    data_lose = data_filt(sad_samples(tt)+wind,:);
    Cal(:,:,tt) = cov(data_lose);
end
Cal = mean(Cal,3);


%% Beamformer
grid = dataprep.leadfield; 

L = grid.leadfield;
W = cell(size(L));

Zstat_mood(1:size(L,2)) = {0};
Tstat_mood(1:size(L,2)) = {0};

for ii = 1:length(L)
    lf = L{ii};
    if ~isempty(lf)
        
        % %           G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        if d(3) < 1
            jj = 2; % The minumum singular value is degenerate
        else
            jj =3;
        end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        w = Cr\lfo / (lfo'/Cr*lfo) ;
        
        %             Qap = w'*Cap*w;
        %             Qan = w'*Can*w;
        
        Qah = w'*Cah*w;
        Qal = w'*Cal*w;
        
        n = w'*noise*w;
        
        W{ii} = w;
        %             Tstat{ii} = (Qa - Qc) ./ (na + nc);
        %             Tstat_rpe{ii} = (Qap - Qan) ./ (Qap + Qan); % normalised version
        Tstat_mood{ii} = (Qah - Qal) ./ (Qah + Qal); % normalised version
        Zstat_mood{ii} = ((Qah + Qal)/2 ) ./ n;
        
    end
    
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
    
end
clc
fprintf('SAM finished\n' )
%% Save result
mri = dataprep.mri;

if ~exist([data_path,'results'],'dir')
    mkdir results
end
cd results

mriname = 'ft_coreg_anat.nii';
if ~exist(mriname,'file')
    ft_write_mri(mriname,mri,'dataformat','nifti');
end


cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';

%     if ~exist(zname,'file')

sourceTstat.avg.pow =  cell2mat(Zstat_mood);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

zname =  sprintf('MoodZ_3s_PseudoT_%d-%dHz_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');


sourceTstat.avg.pow =  cell2mat(Tstat_mood);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

% 3s before mood rating high - low mood
zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');

%%
%     crang = [];
%     sourcePostInt.anatomy = mri.anatomy;
%     cfg = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter = 'pow';
%     cfg.maskparameter = 'pow';
%     cfg.funcolormap  = 'auto';
%     cfg.funcolorlim   = crang;
%     cfg.opacitylim = crang;
%     ft_sourceplot(cfg, sourcePostInt);

%% Template transform
% cfg = [];
% cfg.parameter = 'anatomy';
% [mrinorm] = ft_volumenormalise(cfg, mri);
sourcePostInt.anatomy = mri.anatomy;
cfg = [];
cfg.parameter = 'pow';
[sourcenorm] = ft_volumenormalise(cfg, sourcePostInt);

if ~exist('ft_coreg_anat_norm.nii','file')
    zname =  sprintf('ft_coreg_anat_norm');
    zname_nii = [zname,'.nii'];
    ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');
end

sourcenorm.anatomy = sourcenorm.pow;
zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz_norm_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');

