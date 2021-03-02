clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
data_path = '/data/MBDU/MEG_Booster/data/derivatives/sub-23732/';
cd(data_path)

mri_name = 'anat/23732_anat_visit5_1+orig.HEAD';
% mri_name = 'anat/ortho+orig.HEAD';

data_name = 'meg/23732_gambling_20190515_proc.ds';

mri = ft_read_mri(mri_name);

tagset_shape = mri.hdr.TAGSET_NUM;
tagset_coord = mri.hdr.TAGSET_FLOATS;
tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa
tagset_coord = tagset_coord([2,3,1],:)'; % fiducials have shuffled coordinates
mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin 

fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';
fiducial_coord = abs([257, 255, 0, 0] - fiducial_coord);

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

% figure; 
% subplot(131); imagesc(squeeze(mri.anatomy(238,:,:))); set(gca,'YDir','normal')
% subplot(132); imagesc(squeeze(mri.anatomy(:,131,:))); set(gca,'YDir','normal')
% subplot(133); imagesc(squeeze(mri.anatomy(:,:,100))); set(gca,'YDir','normal')


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

sens = ft_read_sens(data_name);

%% Reading MEG data

freqband = cell(1,5);
freqband{1} = [1 4];
freqband{2} = [4 8];
freqband{3} = [8 13];
freqband{4} = [13 30];
freqband{5} = [40 100];

PseudoT = cell(1,5);

for bb = [4,5]
    
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';
    cfg.bpfreq = freqband{bb};
%     cfg.bpfilttype    = 'firls'; % ???? 
%     cfg.bpfiltord     = 200;
    cfg.bpfilter      = 'yes';
    data = ft_preprocessing(cfg);
    
    data_filt = data.trial{1}';
    
    mu = 4;
    C = cov(data_filt);
    noise = min(svd(C)).*eye(size(C));
    Cr = C + mu*noise;
    
%     figure; subplot(121); imagesc(C);
%     title('Data covariance');
%     subplot(122); imagesc(Cr);
%     title('Regularized data covariance');
    
    % Check power spectra
%     nfft = 2^13;
%     F = fft(data_filt, nfft);
%     P = abs(F(1:nfft/2+1, :));
%     freqs = linspace(0,data.fsample/2, nfft/2+1);
%     figure(1); clf; plot(freqs, P)
    
    %%
    
    cfg                 = [];
    cfg.grad            = sens;
    cfg.headmodel       = vol;
    cfg.reducerank      = 2;
    cfg.channel         = {'MEG'};
    cfg.resolution = 0.5;   % use a 3-D grid with a 0.5 cm resolution
    cfg.sourcemodel.unit       = 'cm';
    cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
    [grid] = ft_prepare_leadfield(cfg);
    
    %% Beamformer
    
    event = ft_read_event(data_name);
    
    button_sample = cell(size(event));
    baseline_sample = cell(size(event));
    
    wind_samples = 1*data.fsample; % 1s window
    
    
    % Try localizing button presses in all frequencies to check it is
    % working correctly!!
    for tt = 1:length(event)
        if strcmp(event(tt).type, 'UPPT001' )
            if event(tt).value == 16 || event(tt).value == 17 % boost win
                button_sample{tt} = event(tt).sample + (0: wind_samples-1);
            
                baseline_sample{tt} = event(tt).sample + (-2*wind_samples: -wind_samples-1);
            end
            
        end
    end
    
    baseline_sample = cell2mat(baseline_sample);
    button_sample = cell2mat(button_sample);
    
    Cc = zeros([size(C),size(baseline_sample,1)]);
    for tt= 1:size(baseline_sample,1)
        data_baseline = data_filt(baseline_sample(tt,:),:);
        Cc(:,:,tt) = cov(data_baseline);
    end
    Cc = mean(Cc(:,:,randperm(size(baseline_sample,1),size(button_sample,1))),3);
    
    
    Ca = zeros([size(C),size(button_sample,1)]);
    for tt= 1:size(button_sample,1)
        data_button = data_filt(button_sample(tt,:),:);
        Ca(:,:,tt) = cov(data_button);
    end
    Ca = mean(Ca,3);
    
    Cna = eye(size(Ca))*min(svd(Ca));
    Cnc = eye(size(Cc))*min(svd(Cc));
    %%
    L = grid.leadfield;
    W = cell(size(L));
    Tstat = cell(size(L));
    
    dipout = cell(size(L));
    
    
    for ii = 1:length(L)
        lf = L{ii};
        if ~isempty(lf)
            
            % fieldtrip method
%             [u, s, v] = svd(real(pinv(lf' / (Cr) *lf)));
%             eta = u(:,1);
%             lfo  = lf * eta;
%             
%             dipout{ii} = eta;
             
%           G O'Neill method, equivalent to ft              
            [v,d] = svd(lf'/Cr*lf);
            d = diag(d);
            if d(3) < 1
                jj = 2; % The minumum singular value is degenerate
            else
                jj =3; 
            end
            lfo = lf*v(:,jj); % Lead field with selected orientation
            
            dipout{ii} = v(:,jj);
            
            w = Cr\lfo / (lfo'/Cr*lfo) ;
            
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
    
    PseudoT{bb} = cell2mat(Tstat);
    
end
%%

bb =4;
cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow = PseudoT{bb};

sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);

crang = [-15 -5];

cfg = [];
cfg.method        = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolorlim   = crang;
cfg.opacitymap = [];
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcePostInt);
set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)))


