clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490'];

s = 5 % [1,3:4]


sub = subn(s,:); % with pixel '24138'artefacts; '24103';    % no pixel '24172'; '24071'
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 1-150 Hz to adjust baseline
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
cfg.viewresult = 'yes';

mri = ft_volumerealign(cfg,mri);

%% Segment MRI
cfg = [];
cfg.output  = 'brain';
segmentmri = ft_volumesegment(cfg,mri);

%% Head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentmri);
sens = ft_read_sens(data_name,'senstype','meg');

%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.resolution      = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit   = 'cm';
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)

if ~exist([data_path,'/lead_fields.mat'],'file')
    [grid] = ft_prepare_leadfield(cfg);
    save([data_path,'/lead_fields'],'grid')
else
    load([data_path,'/lead_fields.mat'])
end
%% Clean Data

% Unsure of best pre-processing for baseline in sub 24071
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
% cfg.demean = 'yes';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [1 150];
data = ft_preprocessing(cfg);


if exist([data_path,'ICA_artifacts.mat'],'file')
    load([data_path,'ICA_artifacts.mat']);
else
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = 'EEG';
    % cfg.demean = 'yes';
    % cfg.bpfilter = 'yes';
    % cfg.bpfreq = [1 150];
    eog = ft_preprocessing(cfg);
    eog = eog.trial{1}(1,:);
    
    
    cfg =[];
    cfg.method = 'pca';
    comp_pca = ft_componentanalysis(cfg, data);
    score = comp_pca.trial{1}';
    compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
    icomp = nnz(compvar95) ;
    clc
    fprintf('%d components for 95perc. of data variance\n',icomp)
    
    if icomp>30
        icomp = 30;
    end
    cfg =[];
    cfg.method = 'fastica';
    cfg.fastica.numOfIC = icomp;
    comp = ft_componentanalysis(cfg, data);
    
    
    figure
    cfg           = [];
    cfg.component = [1:icomp];       % specify the component(s) that should be plotted
    cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)
    
    
    cfg          = [];
    cfg.channel  = [1:5]; % components to be plotted
    cfg.viewmode = 'component';
    cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
    ft_databrowser(cfg, comp)
    
    figure;
    plot(abs(corr(eog',comp.trial{1}')))
    xlabel('ICA component')
    ylabel('Correlation with EOG')
    
    icadel = input('ICA component to eliminate (input as [''01'';''02'']): ');
    
    cfg = [];
    cfg.channel = cell(size(icadel,1),1);
    for ii = 1:size(icadel,1)
        cfg.channel{ii}  = ['fastica0',icadel(ii,:)];
    end
    
    [comps] = ft_selectdata(cfg, comp);
    save([data_path,'/ICA_artifacts'],'comps')
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data_clean    = ft_rejectcomponent(cfg, comps,data);


%% Read events

bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');
for ii = 1:length(bv_names)
    if strncmp(bv_names(ii).name,sub,5)
        bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
bv = readtable(bv_name,opts); % bahavioral data

ratemood_block = bv.blockHappyText_started;
ratemood_block(isnan(ratemood_block))= [] ;


% pix = ft_read_data(data_name,'chanindx','UADC016');
% Outdated way to read light pixel channel
% cfg = [];
% cfg.dataset = data_name;
% cfg.continuous = 'no';
% cfg.channel = 'UADC016';
% cfg.demean = 'yes';
% pix = ft_preprocessing(cfg);
%
% pixm  = cell2mat(pix.trial)';
% pix_white = pixm>2.5;
% pix_sample = find(diff(pix_white)==1);

% read LIGHT marker
ii = 0;
event = ft_read_event(data_name);
pix_sample = zeros(size(event));
for tt = 1:length(event)
    if strcmp(event(tt).type, 'LIGHT' )
        ii = ii +1;
        pix_sample(ii) = event(tt).sample;
    end
end
d = diff(pix_sample);
ii = find(d<(0.25*1200));
pix_sample(ii+1) = []; % correctly identifies triggers and eliminates zeros

if isempty(pix_sample)
    ii = 0;
    event = ft_read_event(data_name);
    for tt = 1:length(event)
        if strcmp(event(tt).type, 'UPPT001' )
            ii = ii +1;
            pix_sample(ii) = event(tt).sample + 24; % Average 24 sample delay
        end
    end
end
f =data.fsample;

bv_answer = bv.fixCross_started;
bv_answer(isnan(bv_answer))=[];

bv_choice = bv.fixCross_2_started;
bv_choice(isnan(bv_choice)) = [];

bv_outcome = bv.fixCross_3_started;
bv_outcome(isnan(bv_outcome)) = [];

bv_mood_block = bv.blockHappyText_started;
bv_mood_block(isnan(bv_mood_block)) = [];

bv_slider_block = bv.blockHappySlider_started;
bv_slider_block(isnan(bv_slider_block)) = [];

bv_mood = bv.happyText_started;
bv_mood(isnan(bv_mood)) = [];

bv_slider = bv.happySlider_started;
bv_slider(isnan(bv_slider)) = [];

bv_rest = bv.endOfBlockText_started;
bv_rest(isnan(bv_rest)) = [];

% Check time difference between pixel and behavioral file
bv_all = sort([bv_answer; bv_choice; bv_outcome; bv_mood; bv_mood_block; bv_slider; bv_slider_block; bv_rest]);

bv_timeoffset = median(bv_all(1:length(pix_sample))-pix_sample/f); % difference in start time between the behavioral and MEG data.
if std(bv_all(1:length(pix_sample))-pix_sample/f) > 20e-3
    warning('MEG and behavioral file temporal mismatch > 20ms')
    figure; histogram(bv_all(1:length(pix_sample))-pix_sample/f);
    xlabel('Time off-set (seconds)');
end


answer_match = [];
choice_match = [];
outcome_match =[];
ITI_match = [];
mood_match = [];
blockmood_match = [];
slider_match = [];
blockslider_match = [];
% Find all choices
nn = 0;
for ii = 1:length(pix_sample)
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        answer_match.pix_index(nn) = ii;
        answer_match.bv_index(nn) = jj;
        answer_match.sample(nn) = pix_sample(ii);
        %        answer_match.RPE(nn) = bv.RPE(jj);
        %        answer_match.winAmount = bv.winAmount(jj);
    end
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_2_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        choice_match.pix_index(nn) = ii;
        choice_match.bv_index(nn) = jj;
        choice_match.sample(nn) = pix_sample(ii);
        if strcmp(bv.choice{jj},'gamble')
            choice_match.gamble(nn) = 1;
        elseif strcmp(bv.choice{jj},'certain')
            choice_match.gamble(nn) = 0;
        end
    end
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_3_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        outcome_match.pix_index(nn) = ii;
        outcome_match.bv_index(nn) = jj;
        outcome_match.sample(nn) = pix_sample(ii);
        if strcmp(bv.outcome{jj},'win')
            outcome_match.win(nn) = 1;
        elseif strcmp(bv.outcome{jj},'lose')
            outcome_match.win(nn) = -1;
        elseif strcmp(bv.outcome{jj},'certain')
            outcome_match.win(nn) = 0;
        end
    end
    
    % Take fixation cross at the end of each trial as baseline
    % lasts 2s at the end of each trial
    sdiff = abs(pix_sample(ii)-2*f - (bv.fixCross_ITI_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        ITI_match.pix_index(nn) = ii;
        ITI_match.bv_index(nn) = jj;
        ITI_match.sample(nn) = pix_sample(ii)-2*f;
        %             if strcmp(bv.outcome{jj},'win')
        %                 ITI_match.win(nn) = 1;
        %             elseif strcmp(bv.outcome{jj},'lose')
        %                 outcome_match.win(nn) = -1;
        %             elseif strcmp(bv.outcome{jj},'certain')
        %                 outcome_match.win(nn) = 0;
        %             end
    end
    
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.happySlider_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        slider_match.pix_index(nn) = ii;
        slider_match.bv_index(nn) = jj;
        slider_match.sample(nn) = pix_sample(ii);
        slider_match.mood(nn) = bv.happySlider_response(jj);
    end
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.blockHappySlider_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        blockslider_match.pix_index(nn) = ii;
        blockslider_match.bv_index(nn) = jj;
        blockslider_match.sample(nn) = pix_sample(ii);
        blockslider_match.mood(nn) = bv.blockHappySlider_response(jj);
    end
    
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.happyText_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        mood_match.pix_index(nn) = ii;
        mood_match.bv_index(nn) = jj;
        mood_match.sample(nn) = pix_sample(ii);
        mood_match.mood(nn) = bv.happySlider_response(jj);
    end
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.blockHappyText_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        blockmood_match.pix_index(nn) = ii;
        blockmood_match.bv_index(nn) = jj;
        blockmood_match.sample(nn) = pix_sample(ii);
        blockmood_match.mood(nn) = bv.blockHappySlider_response(jj);
    end
    
    
end


moodblock = blockmood_match.mood(blockmood_match.mood>0);
mood0 = moodblock(1);

pRPE = outcome_match.win == 1 ;
nRPE = outcome_match.win == -1 ;
pRPE_sample = outcome_match.sample(pRPE);
nRPE_sample = outcome_match.sample(nRPE);

%% Filter data in theta band
freqband = [4,8];

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, freqband,filt_order,'but')';

mu = 4;
C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));
Cr = C + mu*noise;


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
% Cap = mean(Cap,3);
% 
% Can = zeros([size(C),size(lose_samples,2)]);
% for tt= 1:size(lose_samples,2)
%     data_lose = data_filt(lose_samples(tt)+wind,:);
%     Can(:,:,ii) = cov(data_lose);
% end
% Can = mean(Can,3);
% 



%% Mood decision 3s

ind = mood_match.sample>0;
mood_match.mood(~ind)= NaN;

v1 = mood_match.mood(ind);

ind = blockmood_match.sample>0;
blockmood_match.mood(~ind) = NaN;
v0 = blockmood_match.mood(ind);

v = [v1,v0];

highmood = mood_match.mood > (mean(v)+std(v)/2 );
lowmood = mood_match.mood < (mean(v)-std(v)/2 );

win_samples = mood_match.sample(highmood);
lose_samples = mood_match.sample(lowmood);

highmood = blockmood_match.mood > (mean(v)+std(v)/2 );
lowmood = blockmood_match.mood < (mean(v)-std(v)/2 );

win_samples = [win_samples, blockmood_match.sample(highmood)];
lose_samples = [lose_samples, blockmood_match.sample(lowmood)];

wind = 0.25*f:2.75*f; % Exclude edges for sensory and motor overlap
% cuep3_sample = zeros(nnz(choice_match.sample), length(wind));


Cah = zeros([size(C),size(win_samples,2)]);
for tt= 1:size(win_samples,2)
    data_win = data_filt(win_samples(tt)+wind,:);
    Cah(:,:,ii) = cov(data_win);
end
Cah = mean(Cah,3);

Cal = zeros([size(C),size(lose_samples,2)]);
for tt= 1:size(lose_samples,2)
    data_lose = data_filt(lose_samples(tt)+wind,:);
    Cal(:,:,ii) = cov(data_lose);
end
Cal = mean(Cal,3);


%% Beamformer
L = grid.leadfield;
W = cell(size(L));

Tstat_rpe(1:size(L,2)) = {0};
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
        
        Qap = w'*Cap*w;
        Qan = w'*Can*w;
        
        Qah = w'*Cah*w;
        Qal = w'*Cal*w;
        
        W{ii} = w;
        %             Tstat{ii} = (Qa - Qc) ./ (na + nc);
        Tstat_rpe{ii} = (Qap - Qan) ./ (Qap + Qan); % normalised version
        Tstat_mood{ii} = (Qah - Qal) ./ (Qah + Qal); % normalised version
        
    end
    
    if mod(ii,100) == 0
        clc
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
    
end
clc
fprintf('SAM finished\n' )
%% Save result

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

sourceTstat.avg.pow =  cell2mat(Tstat_rpe);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

zname =  sprintf('RPEdiff_0.2s_PseudoT_%d-%dHz',freqband(1),freqband(2));
zname_nii = [zname,'.nii'];
ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');


sourceTstat.avg.pow =  cell2mat(Tstat_mood);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;

% 3s before mood rating high - low mood
zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz',freqband(1),freqband(2));
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
[mrinorm] = ft_volumenormalise(cfg, mri);
[sourcenorm] = ft_volumenormalise(cfg, sourcePostInt);
%% Compare mood
x = mood_match.sample(mood_match.sample>0)';
x1 =  blockmood_match.sample( blockmood_match.sample > 0)';
v = mood_match.mood(mood_match.sample>0)';
v1 = blockmood_match.mood(blockmood_match.sample>0)';
[x,ii] = sort([x ; x1]);
v = [v; v1 ];
v = v(ii);

figure
plot(x/f, v)


F = griddedInterpolant(x,v,'pchip');
xi = x(1):length(data_clean.time{1});
vi = F(xi); % Mood timecourse

hold on
plot(xi/f,vi)