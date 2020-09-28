clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
sub = '24085'; % '23732';
nv ='1'; % '5'
data_path = ['/data/MBDU/MEG_Booster/data/derivatives/sub-',sub,'/'];
cd(data_path)

mri_name = ['anat/',sub,'_anat_visit',nv,'_1+orig.BRIK'];
data_name = ['meg/',sub,'_gambling_20190520_proc.ds'];

mri = ft_read_mri(mri_name);

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
% subplot(131); imagesc(squeeze(mri.anatomy(238,:,:))); set(gca,'YDir','normal')
% subplot(132); imagesc(squeeze(mri.anatomy(:,131,:))); set(gca,'YDir','normal')
% subplot(133); imagesc(squeeze(mri.anatomy(:,:,100))); set(gca,'YDir','normal')
% colormap(grey)

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
cfg.resolution      = 1;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit       = 'cm';
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);


%%

freqband = cell(1,5);
freqband{1} = [1 4];
freqband{2} = [4 8];
freqband{3} = [8 13];
freqband{4} = [13 30];
freqband{5} = [30 60];

PseudoT = cell(1,5);
PseudoZa = cell(1,5);
PseudoZc = cell(1,5);

for bb = 1:5
    
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';
    cfg.bpfreq = freqband{bb};
    
    cfg.bpfilttype    = 'firls'; % fieldtrip reccommends butterworth filter not firls
    cfg.bpfiltord     = 200;
    
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
    
    %% Reading MEG data
    
    
    event = ft_read_event(data_name);
    
    boost_sample = cell(size(event));
    lose_sample = cell(size(event));
    
    wind_samples = 1*data.fsample; % 1s window
    
    
    % Try localizing button presses in all frequencies to check it is
    % working correctly!!
    for tt = 1:length(event)
        if strcmp(event(tt).type, 'UPPT001' )
            if event(tt).value == 70 || event(tt).value == 75 || ...
                    event(tt).value == 60 || event(tt).value == 65% all wins
                boost_sample{tt} = event(tt).sample + (0: wind_samples);
            elseif event(tt).value == 55 || event(tt).value == 50  % all losses
                lose_sample{tt} = event(tt).sample + (0: wind_samples);
            end
            
        end
        
        %     if strcmp(event(tt).type, 'UPPT001' )
        %         if event(tt).value == 70 || event(tt).value == 75 % boost win
        %             boost_sample{tt} = event(tt).sample + (0: wind_samples);
        %         elseif event(tt).value == 55 % lose 25 points
        %             lose_sample{tt} = event(tt).sample + (0: wind_samples);
        %         end
        %
        %     end
        
        
    end
    
    lose_sample = cell2mat(lose_sample);
    boost_sample = cell2mat(boost_sample);
    
    if boost_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
        boost_sample(end,:) = [];
    end
    if lose_sample(end,end)> length(data_filt)
        lose_sample(end,:) = [];
    end
    
    
    %% Beamformer
    k = min(size(lose_sample,1),size(boost_sample,1));
    
    Cc = zeros([size(C),size(lose_sample,1)]);
    for tt= 1:size(lose_sample,1)
        data_lose = data_filt(lose_sample(tt,:),:);
        Cc(:,:,tt) = cov(data_lose);
    end
    Cc = mean(Cc(:,:,randperm(size(lose_sample,1),k)),3);
%     Cc = mean(Cc,3);
    
    Ca = zeros([size(C),size(boost_sample,1)]);
    for tt= 1:size(boost_sample,1)
        data_boost = data_filt(boost_sample(tt,:),:);
        Ca(:,:,tt) = cov(data_boost);
    end
%     Ca = mean(Ca,3);
    Ca = mean(Ca(:,:,randperm(size(boost_sample,1),k)),3);
    
    Cna = eye(size(Ca))*min(svd(Ca));
    Cnc = eye(size(Cc))*min(svd(Cc));
    %%
    L = grid.leadfield;
    W = cell(size(L));
    
    dipout = cell(size(L));
    
    Za(1:size(L,2)) = {0};
    Zc(1:size(L,2)) = {0};
    Tstat(1:size(L,2)) = {0};
    
    parfor ii = 1:length(L)
        lf = L{ii};
        if ~isempty(lf)
            
%             fieldtrip method
%             [u, s, v] = svd(real(pinv(lf' / (Cr) *lf)));
%             eta = u(:,1);
%             lfo  = lf * eta;
%             
%             dipout{ii} = eta;
             
% %           G O'Neill method, equivalent to ft              
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
%             Tstat{ii} = (Qa - Qc) ./ (Qa + Qc); % normalised version
            
            Za{ii} = Qa / na;
            Zc{ii} = Qc/ nc;
            
% 
%             wnorm = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo );
%             Qa = wnorm'*Ca*wnorm;
%             Qc = wnorm'*Cc*wnorm;
            
%             Za{ii} = Qa;
%             Zc{ii} = Qc;
        end
    end
   PseudoT{bb} = cell2mat(Tstat);
   PseudoZa{bb} = cell2mat(Za);
   PseudoZc{bb} = cell2mat(Zc);
end 


%%
cd(data_path)

if ~exist('result','dir')
    mkdir('result')
end

cd result
mriname = 'ft_coreg_anat.nii';
if ~exist(mriname,'file')
    ft_write_mri(mriname,mri,'dataformat','nifti');
end

for bb =1:5
       
    zname = sprintf('PseudoZ_win_1s_%d-%dHz.nii',freqband{bb}(1),freqband{bb}(2));
    
    cfg = [];
    cfg.parameter = 'pow';
    sourceTstat = struct;
    sourceTstat.dim = grid.dim;
    sourceTstat.inside = grid.inside;
    sourceTstat.pos = grid.pos;
    sourceTstat.method = 'average';
    
%     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoZa{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);        
        sourcePostInt.anatomy = sourcePostInt.pow;
              
        ft_write_mri(zname,sourcePostInt,'dataformat','nifti');
%     end
    
    zname = sprintf('PseudoZ_lose_1s_%d-%dHz.nii',freqband{bb}(1),freqband{bb}(2));
    
%     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoZc{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        ft_write_mri(zname,sourcePostInt,'dataformat','nifti');
%     end
    
    tname =  sprintf('PseudoT_win-lose_1s_%d-%dHz.nii',freqband{bb}(1),freqband{bb}(2));
        
    cfg = [];
    cfg.parameter = 'pow';
    sourceTstat = struct;
    sourceTstat.dim = grid.dim;
    sourceTstat.inside = grid.inside;
    sourceTstat.pos = grid.pos;
    sourceTstat.method = 'average';
        
    
%     if ~exist(tname,'file')
        
        sourceTstat.avg.pow =  PseudoT{bb};        
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        ft_write_mri(tname, sourcePostInt,'dataformat','nifti');
%     end
    
end
%% Plots
% close all
% bb =1 ;
% 
% cfg = [];
% cfg.parameter = 'pow';
% sourceTstat = struct;
% sourceTstat.dim = grid.dim;
% sourceTstat.inside = grid.inside;
% sourceTstat.pos = grid.pos;
% sourceTstat.method = 'average';
% sourceTstat.avg.pow =  (PseudoZa{bb} + PseudoZc{bb})/2;
% 
% sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
% 
% crang = [100 350];
% 
% cfg = [];
% cfg.method        = 'ortho';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourcePostInt);
% set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
%     'position',[50 600 845 719])
% 
% 
% cfg = [];
% cfg.parameter = 'pow';
% sourceTstat = struct;
% sourceTstat.dim = grid.dim;
% sourceTstat.inside = grid.inside;
% sourceTstat.pos = grid.pos;
% sourceTstat.method = 'average';
% sourceTstat.avg.pow =  PseudoT{bb};
% 
% sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
% 
% crang = [];
% 
% cfg = [];
% cfg.method        = 'ortho';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourcePostInt);
% set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
%     'position',[1000 600 845 719])

