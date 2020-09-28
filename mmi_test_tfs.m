clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
sub = '24071'; %'24138';
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

mri_name = [sub,'_anat_1+orig.BRIK'];
% mri_name = [sub,'_anat_orient+orig.BRIK'];

data_name = [sub,'MMI_mmi3_proc.ds'];

if ~exist(mri_name,'file')
    unix(['gunzip ',mri_name])
end

mri = ft_read_mri(mri_name,'dataformat','afni_brik');

% mri_name = [sub,'_anat_1+orig.BRIK'];
% mrico = ft_read_mri(mri_name,'dataformat','afni_brik');
tagset_shape = mri.hdr.TAGSET_NUM;
tagset_coord = mri.hdr.TAGSET_FLOATS;
tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa

tagset_p = zeros(1,3);  % Ideal orientation {RL; PA; IS}
% m = eye(4);  m(:,4) = 1;
for ii =1:3
    if strcmp(mri.hdr.Orientation(ii,:),'AP') || strcmp(mri.hdr.Orientation(ii,:),'PA')
        tagset_p(ii) = 2;
%         if strcmp(mri.hdr.Orientation(ii,:),'AP')
%         m(ii,ii) = -1;   m(ii,4) = mri.dim(ii);
%         end
    elseif strcmp(mri.hdr.Orientation(ii,:),'LR') || strcmp(mri.hdr.Orientation(ii,:),'RL')
        tagset_p(ii) = 1; 
%         if strcmp(mri.hdr.Orientation(ii,:),'LR')
%         m(ii,ii) = -1;   m(ii,4) = mri.dim(ii);
%         end
    elseif strcmp(mri.hdr.Orientation(ii,:),'SI') || strcmp(mri.hdr.Orientation(ii,:),'IS')
        tagset_p(ii) = 3;
%          if strcmp(mri.hdr.Orientation(ii,:),'SI')
%             m(ii,ii) = -1;   m(ii,4) = mri.dim(ii);
%         end
    end
end

m = [   -1  0   0   mri.dim(1)
        0   -1  0   mri.dim(2)
        0   0   1   1
        0   0   0   1] ;


tagset_coord = tagset_coord(tagset_p,:)'; % fiducials have shuffled coordinates

mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin 

% fiducial_coord = [tagset_coord,ones(3,1)]*(m/(mri.transform ))';

mri.transform = mri.transform/m;
fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';

% 
cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'yes';

mri = ft_volumerealign(cfg,mri);

% figure; 
% subplot(131); imagesc(squeeze(mri.anatomy(fiducial_coord(1,1),:,:))); set(gca,'YDir','normal')
% subplot(132); imagesc(squeeze(mri.anatomy(:,fiducial_coord(1,2),:))); set(gca,'YDir','normal')
% subplot(133); imagesc(squeeze(mri.anatomy(:,:,fiducial_coord(1,3)))); set(gca,'YDir','normal')
% colormap(gray)
% 
% 
% figure; 
% subplot(131); imagesc(squeeze(mri.anatomy(fiducial_coord(2,1),:,:))); set(gca,'YDir','normal')
% subplot(132); imagesc(squeeze(mri.anatomy(:,fiducial_coord(2,2),:))); set(gca,'YDir','normal')
% subplot(133); imagesc(squeeze(mri.anatomy(:,:,fiducial_coord(2,3)))); set(gca,'YDir','normal')
% colormap(gray)
% 
% figure; 
% subplot(131); imagesc(squeeze(mri.anatomy(fiducial_coord(3,1),:,:))); set(gca,'YDir','normal')
% subplot(132); imagesc(squeeze(mri.anatomy(:,fiducial_coord(3,2),:))); set(gca,'YDir','normal')
% subplot(133); imagesc(squeeze(mri.anatomy(:,:,fiducial_coord(3,3)))); set(gca,'YDir','normal')
% colormap(gray)

%% Segment MRI
cfg = [];
cfg.output  = 'brain';
segmentmri = ft_volumesegment(cfg,mri);

%% Head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentmri);

sens = ft_read_sens(data_name);


%% Calculate lead fields




cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
% cfg.resolution      = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit       = 'cm';
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
cfg.sourcemodel.pos = [ 8.1  -.7   3.9   % right frontal 
                        1.4  3.8   6.5     % left amygdala? 
                        2.3  -4.7  6.4    % right amygdala
                        4.6  0.8   4.9 ];
[tfspos] = ft_prepare_leadfield(cfg);


%%

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
cfg.bpfreq = [1 300];
% cfg.bpfilttype    = 'firls'; % fieldtrip reccommends butterworth filter not firls
% cfg.bpfiltord     = 200;
cfg.bpfilter      = 'yes';
data = ft_preprocessing(cfg);


%%

[coeff,score,latent]=pca( data.trial{1}');

compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
icomp = nnz(compvar95) ;
fprintf('%d components for 95perc. of data variance\n',icomp)

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
cfg.channel  = [1:10]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, comp)

FF = fft(comp.trial{1}');
P = abs(FF(1:length(FF)/2+1,:));
figure; plot(linspace(0,600,length(FF)/2+1), P(:,11))

cfg           = [];
cfg.component = [11 21];
data_clean    = ft_rejectcomponent(cfg, comp,data);

data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [1, 150] )';


%% Prepare covariance
mu = 4;
C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));
Cr = C + mu*noise;

%% Reading MEG data


event = ft_read_event(data_name);

win_sample = cell(size(event));
lose_sample = cell(size(event));
mood_sample = cell(size(event));
gamble_sample = cell(size(event));
certain_sample = cell(size(event));
answer_sample = cell(size(event));
wind_samples = 1*data.fsample; % 1s window


% Try localizing button presses in all frequencies to check it is
% working correctly!!
for tt = 1:length(event)
    if strcmp(event(tt).type, 'ratemood' ) % Slider appear 3s after
        mood_sample{tt} = event(tt).sample + (-wind_samples: wind_samples*3);
        
    elseif strcmp(event(tt).type, 'Answer' )
        answer_sample{tt} = event(tt).sample + (-wind_samples: wind_samples*3);
    end
    
    if  strcmp(event(tt).type, 'win' )
        win_sample{tt} = event(tt).sample + (-wind_samples : 3*wind_samples); % -1-3 s before new choice
    elseif  strcmp(event(tt).type, 'gamble' )
        gamble_sample{tt} = event(tt).sample + (-wind_samples : 4* wind_samples); % 4 seconds before outcome appears
    elseif  strcmp(event(tt).type, 'certain' )
        certain_sample{tt} = event(tt).sample + (-wind_samples : 3*wind_samples);
    elseif  strcmp(event(tt).type, 'lose' )
        lose_sample{tt} = event(tt).sample + (-wind_samples : 3*wind_samples);
    end
    
    
    %     if strcmp(event(tt).type, 'UPPT001' )
    %         if event(tt).value == 70 || event(tt).value == 75 % win win
    %             win_sample{tt} = event(tt).sample + (0: wind_samples);
    %         elseif event(tt).value == 55 % lose 25 points
    %             lose_sample{tt} = event(tt).sample + (0: wind_samples);
    %         end
    %
    %     end
    
    
end

mood_sample = cell2mat(mood_sample);
lose_sample = cell2mat(lose_sample);
win_sample = cell2mat(win_sample);
gamble_sample = cell2mat(gamble_sample);
certain_sample = cell2mat(certain_sample);
answer_sample = cell2mat(answer_sample);

if mood_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
    mood_sample(end,:) = [];
end
if lose_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
    lose_sample(end,:) = [];
end
if win_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
    win_sample(end,:) = [];
end
if gamble_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
    gamble_sample(end,:) = [];
end
if answer_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
    answer_sample(end,:) = [];
end

if certain_sample(end,end)> length(data_filt)
    certain_sample(end,:) = [];
end


%% Beamformer


L = tfspos.leadfield;
W = cell(size(L));

dipout = cell(size(L));

VE = zeros(length(data_filt),length(L));

for ii = 1:length(L)
    lf = L{ii};
    if ~isempty(lf)
        
        % %           G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        if d(3) < d(2)*1e-10
            jj = 2; % The minumum singular value is degenerate
        else
            jj =3;
        end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        dipout{ii} = v(:,jj);
        
        w = Cr\lfo / (lfo'/Cr*lfo) ;
        
        VE(:,ii) = data_filt * w;
    end
end

%% TFS, hilbert methosd

tf_range = [1 3; 2  4; 3 5; 4 6; 5 7; 6 8; 8 12 ; ]

for ff = 1:length(tf_range)
    
    hilbert(VE_filt)
end
%%
VE_lose =  VE(lose_sample,:);
VE_lose = reshape(VE_lose, [size(lose_sample), length(L)] );

vef = struct;
vef.time(1:size(lose_sample,1),1) = {linspace(-1,3,size(lose_sample,2))};
for ii = 1:size(lose_sample,1)
vef.trial{ii,1} = squeeze(VE_lose(ii,:,:))';
end
vef.fsample = 1200;
% vef.sampleinfo = [3,size(lose_sample,2) ];
vef.label = {'frontal';'left temp'; 'right temp'; 'VS'};

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'yes';
cfg.foi         = 1:100;
cfg.tapsmofrq    =  log(cfg.foi)+1;
% cfg.t_ftimwin    = 10./cfg.foi;  % 7 cycles per time window
cfg.t_ftimwin    = 3./cfg.foi;  % 7 cycles per time window
cfg.toi          = -1:0.01:3;

freq_mtmc = ft_freqanalysis(cfg,vef);

figure;
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relative';
% cfg.maskstyle    = 'saturation';
cfg.zlim         = [-4 4];%'maxabs';
cfg.ylim        = [0  60];
cfg.interactive  = 'no';
for ii = 1:4
    cfg.channel      = vef.label{ii};
    subplot(2,4,ii)
    ft_singleplotTFR(cfg, freq_mtmc);
end


VE_win =  VE(win_sample,:);
VE_win = reshape(VE_win, [size(win_sample), length(L)] );

vef = struct;
vef.time(1:size(win_sample,1),1) = {linspace(-1,3,size(win_sample,2))};
for ii = 1:size(win_sample,1)
vef.trial{ii,1} = squeeze(VE_win(ii,:,:))';
end
vef.fsample = 1200;
% vef.sampleinfo = [3,size(lose_sample,2) ];
vef.label = {'frontal';'left temp'; 'right temp'; 'VS'};

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'yes';
cfg.foi         = 1:100;
cfg.tapsmofrq    =  log(cfg.foi)+1;
% cfg.t_ftimwin    = 10./cfg.foi;  % 7 cycles per time window
cfg.t_ftimwin    = 3./cfg.foi;  % 7 cycles per time window
cfg.toi          = -1:0.01:3;

freqw_mtmc = ft_freqanalysis(cfg,vef);


cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relative';
% cfg.maskstyle    = 'saturation';
cfg.zlim         = [-4 4];%'maxabs';
cfg.ylim        = [0  60];
cfg.channel      = 'right temp';
cfg.interactive  = 'no';

for ii = 1:4
    cfg.channel      = vef.label{ii};
    subplot(2,4,ii+4)
    ft_singleplotTFR(cfg, freqw_mtmc);
end



% 
% cfg              = [];
% cfg.output       = 'pow';
% cfg.channel      = 'all';
% cfg.method       = 'tfr';
% cfg.keeptrials   = 'yes';
% cfg.foilim       = [1 100];
% cfg.width         = 2;  % 7 cycles per time window
% cfg.toi          = -1:0.05:3;
% 
% freq_wavelet = ft_freqanalysis(cfg,vef);
% 
% 
% cfg              = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'relative';
% cfg.maskstyle    = 'saturation';
% cfg.zlim         = 'maxabs';
% cfg.channel      = 'all';
% cfg.interactive  = 'no';
% figure
% ft_singleplotTFR(cfg, freq_wavelet);

% data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [4, 8] )';
% figure; plot(linspace(0,600,length(data_filt)/2+1), abs(FT(1:length(data_filt)/2+1,1)))


%%

VE_lose =  VE(lose_sample,:);
VE_lose = reshape(VE_lose, [size(lose_sample), length(L)] );

vef = struct;
vef.time(1:size(lose_sample,1),1) = {linspace(-1,3,size(lose_sample,2))};
for ii = 1:size(lose_sample,1)
vef.trial{ii,1} = squeeze(VE_lose(ii,:,:))';
end
vef.fsample = 1200;
% vef.sampleinfo = [3,size(lose_sample,2) ];
vef.label = {'frontal';'left temp'; 'right temp'; 'VS'};



VE_win =  VE(win_sample,:);
VE_win = reshape(VE_win, [size(win_sample), length(L)] );

vef = struct;
vef.time(1:size(win_sample,1),1) = {linspace(-1,3,size(win_sample,2))};
for ii = 1:size(win_sample,1)
vef.trial{ii,1} = squeeze(VE_win(ii,:,:))';
end
vef.fsample = 1200;
% vef.sampleinfo = [3,size(lose_sample,2) ];
vef.label = {'frontal';'left temp'; 'right temp'; 'VS'};


