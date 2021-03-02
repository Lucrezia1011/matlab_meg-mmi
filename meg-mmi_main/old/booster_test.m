clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
data_path = '/data/MBDU/MEG_Booster/data/derivatives/sub-23732/';
cd(data_path)

mri_name = 'anat/23732_anat_visit5_1+orig.BRIK';
data_name = 'meg/23732_gambling_20190515_proc.ds';

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas    = [257-19 131 87]; %position of nasion
cfg.fiducial.lpa    = [257-109 93 13]; %position of LPA
cfg.fiducial.rpa    = [257-109 90 162]; %position of RPA
cfg.viewresult = 'yes';


mri = ft_read_mri(mri_name);
mri = ft_volumerealign(cfg,mri);

% mri_test = ft_read_mri(mri_name);
% 
% imagesc(squeeze(mri.anatomy(:,90,:))); 
% set(gca,'YDir','normal')

%% Segment MRI
cfg = [];
cfg.output  = 'brain';
segmentmri = ft_volumesegment(cfg,mri);

%% Head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentmri);

% Visualize result
volcm = ft_convert_units(vol,'cm');
sens = ft_read_sens(data_name);

figure;
ft_plot_sens(sens);
hold on
ft_plot_headmodel(volcm,'edgealpha',0.5,'facealpha',0.5, 'edgecolor','blue');
set(gcf,'color','w')

% Create single sphere model for comparison
% cfg = [];
% cfg.method = 'singlesphere';
% vols = ft_prepare_headmodel(cfg, segmentmri);
% 
% hold on 
% vols = ft_convert_units(vols,'cm');
% ft_plot_headmodel(vols,'edgealpha',0.5,'facealpha',0.5, 'edgecolor','red');
%% Reading MEG data

% hdr = ft_read_header(data_name);
% data = ft_read_data(data_name);

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
cfg.bpfreq = [13 30];
% cfg.bpfilttype    = 'firls';
% cfg.bpfiltord     = 200;
cfg.bpfilter      = 'yes';
data = ft_preprocessing(cfg);

data_filt = data.trial{1}';

mu = 4;
C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));
Cr = C + mu*noise;

figure; subplot(121); imagesc(C);
title('Data covariance');
subplot(122); imagesc(Cr); 
title('Regularized data covariance');

% Check power spectra
nfft = 2^13;
F = fft(data_filt, nfft);
P = abs(F(1:nfft/2+1, :));
freqs = linspace(0,data.fsample/2, nfft/2+1);
figure(1); clf; plot(freqs, P)



%% Redifing the trials
cfg =[];
cfg.dataset = data_name;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'UPPT001';
cfg.trialdef.eventvalue     = [70 75]; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = 2; % in11 seconds
cfg.trialdef.poststim       = 1; % in seconds

cfg = ft_definetrial(cfg);

cfg.channel = 'MEG';
cfg.demean = 'yes';
dataFIC = ft_preprocessing(cfg);

cfg = [];
cfg.toilim = [-1 0];
dataPre = ft_redefinetrial(cfg, dataFIC);

cfg.toilim = [1 2];
dataPost = ft_redefinetrial(cfg, dataFIC);
%% Cross spectral density matrix

% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.output    = 'powandcsd';
% cfg.tapsmofrq = 4;
% cfg.foilim    = [2 10];
% freqPre = ft_freqanalysis(cfg, dataPre);
% 
% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.output    = 'powandcsd';
% cfg.tapsmofrq = 4;
% cfg.foilim    = [2 10];
% freqPost = ft_freqanalysis(cfg, dataPost);
% 
% figure;
% plot(freqPost.freq, freqPost.powspctrm);

%%

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.sourcemodel.unit       = 'cm';
cfg.normalize = 'yes'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);

%% Beamformer

event = ft_read_event(data_name);

boost_sample = cell(size(event));
lose_sample = cell(size(event));

wind_samples = 1*data.fsample; % 1s window

for tt = 1:length(event)
   if strcmp(event(tt).type, 'UPPT001' )
       if event(tt).value == 70 || event(tt).value == 75 % boost win
           boost_sample{tt} = event(tt).sample + (0: wind_samples); 
       elseif event(tt).value == 55 % lose 25 points
           lose_sample{tt} = event(tt).sample + (0: wind_samples); 
       end
           
   end
end

lose_sample = cell2mat(lose_sample);
boost_sample = cell2mat(boost_sample);

Cc = zeros([size(C),size(lose_sample,1)]);
for tt= 1:size(lose_sample,1)
    data_lose = data_filt(lose_sample(tt,:),:);
    Cc(:,:,tt) = cov(data_lose);
end
Cc = mean(Cc,3);


Ca = zeros([size(C),size(boost_sample,1)]);
for tt= 1:size(boost_sample,1)
    data_boost = data_filt(boost_sample(tt,:),:);
    Ca(:,:,tt) = cov(data_boost);
end
Ca = mean(Ca,3);

Cna = eye(size(Ca))*min(svd(Ca));
Cnc = eye(size(Cc))*min(svd(Cc));
%%
L = grid.leadfield;
W = cell(size(L));
Tstat = cell(size(L));
for ii = 1:length(L)
    l = L{ii};
    if ~isempty(l)
        s = svd(l'/Cr*l);
        [ss, jj] = sort(s,'ascend'); % sort from smallest to largest
        l = l(:,jj(2)); %Take second smallest simgular value as leadfileds have been already rank reduced?
        w = Cr\l / (l'/Cr*l) ;
        
        Qa = w'*Ca*w;
        Qc = w'*Cc*w;
        
        na = w'*Cna*w;
        nc = w'*Cnc*w;
        
        W{ii} = w;
        Tstat{ii} = (Qa - Qc) ./ (na + nc);
    else
        Tstat{ii} = 0;
    end
end


cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow = cell2mat(Tstat);

sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = [3 8];
cfg.opacitylim = [3 8];
ft_sourceplot(cfg, sourcePostInt);


% sourcePostIntflip = sourcePostInt;
% sourcePostIntflip.pow = -sourcePostIntflip.pow ;
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = [3 8];
% cfg.opacitylim = [3 8];
% ft_sourceplot(cfg, sourcePostIntflip);


%%
cfg              = [];
cfg.method       = 'dics';
cfg.frequency    = 6;
cfg.sourcemodel  = grid;
cfg.headmodel    = vol;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = 0;
sourcePost = ft_sourceanalysis(cfg, freqPost);

cfg = [];
cfg.parameter = 'pow';
sourcePostInt  = ft_sourceinterpolate(cfg, sourcePost , mri);


% Different normalization method
% sourceNAI = sourcePost_nocon;
% sourceNAI.avg.pow = sourcePost_nocon.avg.pow ./ sourcePost_nocon.avg.noise;

% cfg = [];
% cfg.downsample = 2;
% cfg.parameter = 'pow';
% sourceNAIInt = ft_sourceinterpolate(cfg, sourceNAI , mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [5.0 9];
cfg.opacitylim    = [5.0 9];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, sourcePostInt);
