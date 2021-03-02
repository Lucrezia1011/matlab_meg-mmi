function sourcenorm = mmi_mood_beamformer_v2(sub,freqband,mu,dataprep)


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


%% Exclude ERF windows

samples_erf = [];

ind = answer_match.sample>0;
samples_erf = cat(2,samples_erf, answer_match.sample(ind));

ind = choice_match.sample>0;
samples_erf = cat(2,samples_erf, choice_match.sample(ind));

ind = outcome_match.sample>0 & outcome_match.win~=0; % certain choice does not change stimulus
samples_erf = cat(2,samples_erf, outcome_match.sample(ind));

ind = mood_match.sample>0;
samples_erf = cat(2,samples_erf, mood_match.sample(ind));

ind = blockmood_match.sample>0;
samples_erf = cat(2,samples_erf, blockmood_match.sample(ind));

ind = blockslider_match.sample>0;
samples_erf = cat(2,samples_erf, blockslider_match.sample(ind));

ind = slider_match.sample>0;
samples_erf = cat(2,samples_erf, slider_match.sample(ind));

% Add rest
% ind = rest_match.sample>0;
% samples_erf = cat(2,samples_erf, rest_match.sample(ind));

ind = ITI_match.sample>0;
samples_erf = cat(2,samples_erf, ITI_match.sample(ind));


cfg = [];
cfg.dataset = dataprep.dataname;
cfg.continuous = 'yes';
cfg.channel = {'UADC005';'UADC006';'UADC007'}; % Time channel!
buttons = ft_preprocessing(cfg);

buttonsd = diff(buttons.trial{1}');
buttonpress = buttonsd>1.5;
[samples,~] = find(buttonpress);
samples = sort(samples);

f = dataprep.data.fsample;
samples(find(diff(samples)<0.05*f)+1) = [];

% Eliminate -50ms to 200ms around each task change and buttons press

samples_erf = cat(2,samples_erf, samples');

wind = -0.05*f:0.45*f;
w = wind(end)-wind(1);
samples_erf_wind = zeros(length(samples_erf),length(wind));
for tt = 1:length(samples_erf)
    samples_erf_wind(tt,:) = samples_erf(tt)+wind;
end 

samples_erf_wind = sort(samples_erf_wind(:));
samples_erf_wind = unique(samples_erf_wind);
d = diff(samples_erf_wind);
ind = d>1 & d<w; % eliminate gaps of <200ms

samples_extra1 = samples_erf_wind(ind);
ind = [logical(0);ind];
samples_extra2 = samples_erf_wind(ind);
samples_extra = [];
for tt=1:length(samples_extra1)
    samples_extra = cat(2, samples_extra, samples_extra1(tt)+1:samples_extra2(tt)-1);
end
samples_erf_wind = sort(cat(1,samples_erf_wind,samples_extra'));

samples_erf_wind(samples_erf_wind>dataprep.data.sampleinfo(2)) = [];

%% Filter data in theta band
% freqband = [4,8];
data_clean = dataprep.data;
f = data_clean.fsample;
filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1},f, [2 40],filt_order,'but')';

data_all = data_filt;
data_all(samples_erf_wind,:)=[];
C = cov(data_all);
noise = min(svd(C)).*eye(size(C)); % Could estimate from noise recording?

% Cr = C + mu*noise;
Cr = C + mu*max(svd(C)).*eye(size(C));


%% Mood decision 3s

% ind = mood_match.sample>0;
% mood_match.mood(~ind)= NaN;
% 
% v1 = mood_match.mood(ind);
% samples1 = mood_match.sample(ind);
% 
% ind = blockmood_match.sample>0;
% blockmood_match.mood(~ind) = NaN;
% v0 = blockmood_match.mood(ind);
% samples0 = blockmood_match.sample(ind);
% 
% v = [v1,v0];
% [s,ind] = sort(v);

x = mood_match.sample(mood_match.sample>0)';
x1 =  blockmood_match.sample( blockmood_match.sample > 0)';
v = mood_match.mood(mood_match.sample>0)';
v1 = blockmood_match.mood(blockmood_match.sample>0)';
[x,ii] = sort([x ; x1]);
v = [v; v1 ];
v = v(ii);
F = griddedInterpolant(x,v,'pchip');
xi = x(1):length(data_clean.time{1});
vi = F(xi); % Mood timecourse



[vs, ind] = sort(vi);
q = floor(length(vi)/4);
sad_samples = ind(1:q);
happy_samples = ind(q*3:q*4);

sad_samples  = sort(sad_samples);
happy_samples = sort(happy_samples);

% figure
% plot(xi/f,vi)
% hold on
% plot(xi(sad_samples)/f,vi(sad_samples),'b.','linewidth',2)
% plot(xi(happy_samples)/f,vi(happy_samples),'r.','linewidth',2)

[~,ind,~] =intersect(sad_samples, samples_erf_wind);
sad_samples(ind)  =[];

[~,ind,~] =intersect(happy_samples, samples_erf_wind);
happy_samples(ind)  =[];

% figure
% plot(xi/f,vi)
% hold on
% plot(xi(sad_samples)/f,vi(sad_samples),'b.','linewidth',2)
% plot(xi(happy_samples)/f,vi(happy_samples),'r.','linewidth',2)

data_filt = ft_preproc_bandpassfilter(data_clean.trial{1},f, freqband,filt_order,'but')';

data_win = data_filt(happy_samples,:);
Cah = cov(data_win);

data_lose = data_filt(sad_samples,:);
Cal = cov(data_lose);
%
% samples = [samples1, samples0];
% samples = samples(ind);
%
% hs = s > median(s);
% ls = s < median(s);
%
% happy_samples = samples(hs);
% sad_samples = samples(ls);
%
% if nnz(hs)> nnz(ls)
%     happy_samples(1) = [];
% elseif nnz(hs)< nnz(ls)
%     sad_samples(end) = [];
% end
%
% wind = 0*f:3*f; % Exclude edges for sensory and motor overlap
% cuep3_sample = zeros(nnz(choice_match.sample), length(wind));


% Cah = zeros([size(C),size(happy_samples,2)]);
% for tt= 1:size(happy_samples,2)
%     data_win = data_filt(happy_samples(tt)+wind,:);
%     Cah(:,:,tt) = cov(data_win);
% end
% Cah = mean(Cah,3);
%
% Cal = zeros([size(C),size(sad_samples,2)]);
% for tt= 1:size(sad_samples,2)
%     data_lose = data_filt(sad_samples(tt)+wind,:);
%     Cal(:,:,tt) = cov(data_lose);
% end
% Cal = mean(Cal,3);


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
cd([data_path,'results'])

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

zname =  sprintf('MoodZ_quartile_PseudoT_%d-%dHz_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');


sourceTstat.avg.pow =  cell2mat(Tstat_mood);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

% 3s before mood rating high - low mood
zname =  sprintf('Mooddiff_quartile_PseudoT_%d-%dHz_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');

%%
    crang = [];
    sourcePostInt.anatomy = mri.anatomy;
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    ft_sourceplot(cfg, sourcePostInt);

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
zname =  sprintf('Mooddiff_quartile_PseudoT_%d-%dHz_norm_mu%.0s',freqband(1),freqband(2),mu);
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcenorm,'dataformat','nifti');

