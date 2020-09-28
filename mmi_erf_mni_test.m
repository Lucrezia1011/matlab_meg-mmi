clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

% Can do the same for ERF!
sn = 9; %[1:9,11,12,14:16]% co-register 2 then run these! %[1,3,4,6,7,8,9,11,14,15,16]   %Subjects showing enough variation in mood

sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline

name_list = dir;
data_names = cell(1);
jj = 0;
for ii = 1:length(name_list)
    if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds')
        jj = jj+1;
        data_names{jj} = name_list(ii).name;
    end
end

runs = 1; %:length(data_names)

%%
data_name = data_names{runs};

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

% if exist([data_path,'results/ft_coreg_anat.nii'],'file')
%     mri = ft_read_mri([data_path,'results/ft_coreg_anat.nii']);
%     mri.coordsys = 'ctf';
% else
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

%% Calculate lead fields

if ~exist([processing_folder,'/leadfields_5mm.mat'],'file')
    cfg                 = [];
    cfg.grad            = sens;
    cfg.headmodel       = vol;
    cfg.reducerank      = 2;
    cfg.channel         = {'MEG'};
    cfg.resolution      = 0.005;                   
    cfg.sourcemodel.unit   = 'm';
    cfg.siunits         = true;
    cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
    [grid] = ft_prepare_leadfield(cfg);
    save([processing_folder,'/leadfields_5mm.mat'],'grid')
else
    load([processing_folder,'/leadfields_5mm.mat']);
end

%% Clean data

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
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data_clean    = ft_rejectcomponent(cfg, comps,data);
end

data = data_clean;
clear data_clean

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [1 150],filt_order,'but');
% data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 );

data.trial{1} = data_filt;
clear data_filt

%% Read events

bv_match = match_triggers_fc(data_name);

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;
buttonpress = bv_match.buttonpress;
%
% pRPE = outcome_match.win == 1 ;
% nRPE = outcome_match.win == -1 ;
% pRPE_sample = outcome_match.sample(pRPE);
% nRPE_sample = outcome_match.sample(nRPE);


%% Beta desynch beamformer

data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [13 30],filt_order,'but');
data_beta = data;
data_beta.trial{1} = data_filt;
clear data_filt


icacomps = length(data.cfg.component);
C = cov(data_beta.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
Cr = C + 4*noiseC; % 10*min for 23911

data_erd = define_trials(buttonpress,data_beta,bv_match.time,[0.1 0.7]);
data_ers = define_trials(buttonpress,data_beta,bv_match.time,[2 2.6]);

% data_act = define_trials(outcome_match.sample(outcome_match.win==1),data_beta,bv_match.time,[0 1]);
% data_con = define_trials(outcome_match.sample(outcome_match.win==-1),data_beta,bv_match.time,[0 1]);

ntrials = min([length(data_act.trial),length(data_con.trial)]);

Ca = zeros(length(data_act.label),length(data_act.label),ntrials);
Cc = Ca;
for ii =1:ntrials
    Ca(:,:,ii) = cov(data_act.trial{ii}');
    Cc(:,:,ii) = cov(data_con.trial{ii}');
end



Ca = mean(Ca,3);
Cc = mean(Cc,3);

L = grid.leadfield;
T(1:size(L,2)) = {0};

parfor ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        %         if d(3) < 1
        %             jj = 2; % The minumum singular value is degenerate
        %         else
        %             jj =3;
        %         end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        w = Cr\lfo / (lfo'/Cr*lfo) ;  
        Qa = w'*Ca*w;
        Qc = w'*Cc*w;
        T{ii} = (Qa - Qc) / (Qa+Qc);
    end
end


cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';

sourceTstat.avg.pow =  cell2mat(T);
sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
sourcePostInt.anatomy = sourcePostInt.pow;
zname =  sprintf('Buttonpress');
% zname =  sprintf('Gamble_outcome_win-lose');
zname_nii = [zname,'.nii'];
writebrik(zname,sourcePostInt)
%% Beamfomer


cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.sourcemodel.pos = [31.133 -54.858 36.591 % right temporal
                       -14.867 0.142  73.591 %PCC? 
                       71.133  -1.858 40.591 %ACC
                       32.133  -45.858 98.591
                       ]*1e-3; % CTF use ALS coordinate order                 
cfg.sourcemodel.unit   = 'm';
cfg.siunits         = true;
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[locs] = ft_prepare_leadfield(cfg);



% cfg = [];
% cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
% datadown = ft_resampledata(cfg, data);

icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

Cr = C + 10*noiseC; % need to normalise because of ICA
% Cr = C + 0.01*eye(nchans)*E(1);

L = locs.leadfield;

% VEp(1:size(L,2)) = {0};
% VEn(1:size(L,2)) = {0};

VE = cell(1,size(L,2));
W(1:size(L,2)) = {0};


for ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        %         if d(3) < 1
        %             jj = 2; % The minumum singular value is degenerate
        %         else
        %             jj =3;
        %         end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        w = Cr\lfo / (lfo'/Cr*lfo) ;       
        %         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
        %         Better Hall's or normalized weights
        wnorm = w/sqrt(w'*noiseC*w);
        VE{ii} = w'*data.trial{1};
        %         w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        W{ii} = w;
    end
   
end
clc
fprintf('Beamformer finished\n' )
%%
dataerf = data;
dataerf.trial{1} = cell2mat(VE');
% dataerf = define_trials(buttonpress,dataerf,bv_match.time,[-4 7]);
dataerf = define_trials(outcome_match.sample(outcome_match.sample~=0),dataerf,bv_match.time,[-3 4]);
dataerf.label = dataerf.label(1:size(dataerf.trial{1},1));

time =dataerf.time{1};
cfg = [];
cfg.method     = 'wavelet';
cfg.taper      = 'hanning';
cfg.pad        = 'nextpow2';
cfg.toi        = time;
cfg.foi        = 1:40;
% cfg.tapsmpfreq = linspace(1,20,120);
cfg.t_ftimwin  = 3./cfg.foi;
cfg.keeptrials = 'yes';
freq_segmented = ft_freqanalysis(cfg, dataerf);

cfg              = [];
cfg.baseline     = [-3 -1];
cfg.baselinetype = 'relchange';
cfg.maskstyle    = 'saturation';
cfg.zlim         = []*1e-17;
cfg.zlim         = [-0.5 2];
cfg.channel      = dataerf.label{3};
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, freq_segmented);

%%

highpass_TFS = [2 4 6 8 10 13 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130]; %splits all data 
lowpass_TFS = [6 8 10 13 15 20 25 30 35 40 45 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150 ];

highpass_TFS = [1 2 4 6 8 10 13 15 20 25 30 35 40 ]; %splits all data 
lowpass_TFS = [4 6 8 10 13 15 20 25 30 35 40 45 60 ];


windtime = [-3 3];
basetime = [0 2];

windsamples= outcome_match.sample(outcome_match.win==1);
windsamples= buttonpress;

basesamples  =ITI_match.sample;

time = linspace(windtime(1),windtime(2),f*diff(windtime));
ve = cell2mat(VE');
P = zeros(length(highpass_TFS),length(time),size(ve,1));

b = zeros(length(highpass_TFS),size(ve,1));

for freq = 1:length(highpass_TFS)
    data_filt = ft_preproc_bandpassfilter(ve, data.fsample, [highpass_TFS(freq) lowpass_TFS(freq)],200,'firls');
    data_filt = abs(hilbert(data_filt'));
    
    dataerf = data;
    dataerf.trial{1} = data_filt';
    % dataerf = define_trials(buttonpress,dataerf,bv_match.time,[-4 7]);
    datatrial = define_trials(windsamples,dataerf,bv_match.time,windtime);
    baseline = define_trials(basesamples,dataerf,bv_match.time,basetime); 
     
    a = cell2mat(datatrial.trial);
    a = reshape(a,[size(ve,1),length(time),size(datatrial.sampleinfo,1)]);
    
    P(freq,:,:) = mean(a,3)';
    
    a = cell2mat(baseline.trial);
    a = reshape(a,[size(ve,1),length(baseline.time{1}),size(baseline.sampleinfo,1)]);
    b(freq,:) = mean(mean(a,2),3);
    
end

for vox = 4 %[7,9]
figure;
pcolor(time, (highpass_TFS+lowpass_TFS)/2, P(:,:,vox)-b(:,vox));
shading interp; colorbar;
caxis([-2 2]*1e-9);  ylim([2.5 50])
end
%%
pp = mean(freq_segmented.powspctrm(:,3,:,:),1);
pp = squeeze(pp);
ppbase = zeros(size(pp,1),1);
for ii = 1:size(pp,1)
    ppf = pp(ii,:);
    ppbase(ii) = mean(ppf(~isnan(ppf))); 
end
figure; pcolor(pp - ppbase);
shading interp; colorbar



P = squeeze(freq_segmented.powspctrm);
figure; pcolor(dataerf.time{1},freq_segmented.freq,P')
shading interp

begsample = dataerf.sampleinfo(:,1);
endsample = dataerf.sampleinfo(:,2);
time = ((begsample+endsample)/2) / dataerf.fsample;



ve = cell2mat(dataerf.trial');
plot(dataerf.time{1},mean(ve,1));

Rsamples = outcome_match.sample(outcome_match.win~=0);
dataerf = data;
dataerf.trial{1} = VE{1};
dataerf = define_trials(Rsamples,dataerf,bv_match.time,[-1 2]);
ve = cell2mat(dataerf.trial');
plot(dataerf.time{1},mean(ve,1));



ntrials = length(dataerf.time);
nsamples = length(dataerf.time{1});

dataerfc = cell2mat(dataerf.trial);
trialsRPE = outcome_match.RPE(outcome_match.win~=0);
Rwin = trialsRPE>0;
Rgamble = abs(trialsRPE)>5; % big gambles!
Rsafe = abs(trialsRPE)<3; % safe gambles



dataerf = define_trials(buttonpress,data,bv_match.time,[-1 2]);
cfg = [];
cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
dataerf = ft_resampledata(cfg, dataerf);

ntrials = length(dataerf.time);
nsamples = length(dataerf.time{1});
VEpress = cell(1,size(L,2));
dataerfc = cell2mat(dataerf.trial);
for ii = 1:length(L)
    if ~isempty(L{ii})
        w = W{ii};
        
        dataloc = w'*dataerfc;
        %         dataloc = reshape(dataloc, [nsamples, ntrials]);
        
        VEpress{ii} = dataloc;
        
    end
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
end

clear dataerfc


%%

VE = cell2mat(VE');

%%

VEerf = dataerf;

VEerf.label = cell(size(VE,1),1);
for ii = 1:size(VE,1) % create unique channel names
    VEerf.label{ii} = ['VE',num2str(ii)];
end
for ii = 1:ntrials % create unique channel names
    VEerf.trial{ii} = VE(:,((ii-1)*nsamples +1):(ii*nsamples));
end
% VEica = VEwhite;
% VEica.trial{1} = VEall';
cfg =[];
cfg.method = 'icasso';
cfg.icasso.method = 'fastica';
% cfg.icasso.lastEig = 30;
cfg.numcomponent = 25;
% cfg.fastica.numOfIC = 25;
cfg.channel = 'all';
cfg.icasso.Niter = 10;
comp = ft_componentanalysis(cfg, VEerf);
