clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
sub = '24103'; % with pixel '24138'; '24103';    % no pixel '24172'; '24071'
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


cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcePostInt);
set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
    'position',[1000 600 845 719])


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
cfg.resolution      = 0.5;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit   = 'cm';
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



cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
cfg.bpfreq = [1 300];
% cfg.bpfilttype    = 'firls'; % fieldtrip reccommends butterworth filter not firls
% cfg.bpfiltord     = 200;
cfg.bpfilter      = 'yes';
datab = ft_preprocessing(cfg);


% cfg = [];
% cfg.dataset = data_name;
% % cfg.continuous = 'yes';
% cfg.channel = 'MEG';
% cfg.bpfreq = [1 300];
% % cfg.bpfilttype    = 'firls'; % fieldtrip reccommends butterworth filter not firls
% % cfg.bpfiltord     = 200;
% cfg.bpfilter      = 'yes';
% cfg.demean  =       'yes';
% datab = ft_preprocessing(cfg);

%%
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);


[coeff,score,latent]=pca( data.trial{1}');

compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
icomp = nnz(compvar95) ;
fprintf('%d components for 95perc. of data variance\n',icomp)

cfg =[];
cfg.method = 'fastica';
cfg.fastica.numOfIC = icomp;
comp = ft_componentanalysis(cfg, data);
% 
% cfg =[];
% cfg.method = 'runica';
% comp = ft_componentanalysis(cfg, data);
% 
% cfg =[];
% cfg.method = 'runica';
% bcomp = ft_componentanalysis(cfg, datab);

figure
cfg           = [];
cfg.component = [1:27];       % specify the component(s) that should be plotted
cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)


cfg          = [];
cfg.channel  = [13,14,17,18,26]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, comp)

FF = fft(comp.trial{1}');
P = abs(FF(1:length(FF)/2+1,:));
figure; plot(linspace(0,600,length(FF)/2+1), P(:,11))

cfg           = [];
cfg.component = [17];
data_clean    = ft_rejectcomponent(cfg, comp,data);

% a = A(:,cii)*icasig(cii,:);

%%

bb = 2;

data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [4, 8] )';
% figure; plot(linspace(0,600,length(data_filt)/2+1), abs(FT(1:length(data_filt)/2+1,1)))

%
%     cfg = [];
%     cfg.dataset = data_name;
%     cfg.continuous = 'yes';
%     cfg.channel = 'MEG';
%     cfg.bpfreq = freqband{bb};
%
%     cfg.bpfilttype    = 'firls'; % fieldtrip reccommends butterworth filter not firls
%     cfg.bpfiltord     = 200;
%
%     cfg.bpfilter      = 'yes';
%     data = ft_preprocessing(cfg);
%
%     data_filt = data.trial{1}';

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

f =1200;

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

event = ft_read_event(data_name);

pix = ft_read_data(data_name,'chanindx','UADC016');
pix  = squeeze(pix);
pixm = pix-median(pix);
pixm = reshape(pixm,[size(pixm,1)*size(pixm,2),1]);
pix_white = pixm>3;
pix_sample = find(diff(pix_white)==1);


if isempty(pix_sample)  % no pixel
    win_sample = cell(size(event));
    lose_sample = cell(size(event));
    ratemood_sample = cell(size(event));
    slider_sample = cell(size(event));
    gamble_sample = cell(size(event));
    certain_sample = cell(size(event));
    answer_sample = cell(size(event));
    rest_sample = cell(size(event));
    
    ii = 0;
    event_pix = zeros(size(pix_sample));
    event_pix_value =  zeros(size(pix_sample));
    for tt = 1:length(event)
    
        if strcmp(event(tt).type, 'UPPT001' ) % Slider appear 3s after
            ii = ii+1;
            event_pix(ii) = event(tt).sample;
            event_pix_value(ii) = event(tt).value;
    
            switch event(tt).value
                case 1
                    answer_sample{tt} = event(tt).sample;
                case 2
                    certain_sample{tt} = event(tt).sample;
                    ii = ii+1; % pixel appears when gamble outcome would be
                case 4
                    gamble_sample{tt} = event(tt).sample;
                case 8
                    win_sample{tt} = event(tt).sample;
                case 16
                    lose_sample{tt} = event(tt).sample;
                case 32
                    ratemood_sample{tt} = event(tt).sample;
                case 64
                    slider_sample{tt} = event(tt).sample;
                case 128
                    rest_sample{tt} = event(tt).sample;
            end
        end
    
    end
    
    ratemood_sample = cell2mat(ratemood_sample);
    lose_sample = cell2mat(lose_sample);
    win_sample = cell2mat(win_sample);
    gamble_sample = cell2mat(gamble_sample);
    certain_sample = cell2mat(certain_sample);
    answer_sample = cell2mat(answer_sample);
    rest_sample = cell2mat(rest_sample);
    slider_sample =cell2mat(slider_sample);
    
    
else
    
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
        
    end
    
    
    % Reaction time
    % bv.choiceKey_rt
    %
    % % RPE
    % bv.RPE
    %
    % bv.winAmount
    % bv.loseAmount
    % bv.certainAmount
    % bv.choice
    % bv.outcome
    % bv.outcomeAmount
    
    mm = mean(bv.outcomeAmount(outcome_match.sample>0));
    ss = std(bv.outcomeAmount(outcome_match.sample>0));
    
    bigwin = (bv.outcomeAmount-mm) > ss/2;
    bigloss = (bv.outcomeAmount-mm) < -ss/2;
    
    
    for ii = 1:jj
        win_sample = cat(2,outcome_match.sample(bigwin(1:jj) & abs(outcome_match.win)'),...
            choice_match.sample(bigwin(1:jj) & ~outcome_match.win' ) );
        lose_sample =cat(2,outcome_match.sample(bigloss(1:jj) & abs(outcome_match.win)'),...
            choice_match.sample(bigloss(1:jj) & ~outcome_match.win' ) );
    end
    win_amount = bv.outcomeAmount(bigwin(1:jj));
    lose_amount = bv.outcomeAmount(bigloss(1:jj));
    
    
    pRPE = outcome_match.win == 1 ;
    nRPE = outcome_match.win == -1 ;
    pRPE_sample = outcome_match.sample(pRPE);
    nRPE_sample = outcome_match.sample(nRPE);
   
    % RPE column is the estimated RPE 
    
    
    % the presentation of new options (i.e. answer match) should induce a
    % signal representing the probability of winning, based on current mood
    % and values of options, influencing behaviour (gamble/certain choice). 
    
    % How is mood represented?? Change of electrical potential of a set of
    % neurons? Could it change the baseline power in specific bands.
    
end

    %%
    
%     event = ft_read_event(data_name);
%     
%     win_sample = cell(size(event));
%     lose_sample = cell(size(event));
%     ratemood_sample = cell(size(event));
%     gamble_sample = cell(size(event));
%     certain_sample = cell(size(event));
%     answer_sample = cell(size(event));
%     wind_samples = 1*data.fsample; % 1s window
%     % Try localizing button presses in all frequencies to check it is
%     % working correctly!!
%     for tt = 1:length(event)
%         if strcmp(event(tt).type, 'ratemood' ) % Slider appear 3s after
%             ratemood_sample{tt} = event(tt).sample + (-wind_samples: wind_samples*3);
%            
%         elseif strcmp(event(tt).type, 'Answer' ) % Slider appear 3s after
%             answer_sample{tt} = event(tt).sample + (-wind_samples: wind_samples*2);
%         end    
%         
%         if  strcmp(event(tt).type, 'win' )
%             win_sample{tt} = event(tt).sample + (0 : wind_samples); % -1-3 s before new choice
%         elseif  strcmp(event(tt).type, 'gamble' )    
%             gamble_sample{tt} = event(tt).sample + (-wind_samples : 4* wind_samples); % 4 seconds before outcome appears
%         elseif  strcmp(event(tt).type, 'certain' )    
%             certain_sample{tt} = event(tt).sample + (-wind_samples : 3*wind_samples);
%         elseif  strcmp(event(tt).type, 'lose' )   
%             lose_sample{tt} = event(tt).sample + (0 : wind_samples);
%         end
%             
%         
%         %     if strcmp(event(tt).type, 'UPPT001' )
%         %         if event(tt).value == 70 || event(tt).value == 75 % win win
%         %             win_sample{tt} = event(tt).sample + (0: wind_samples);
%         %         elseif event(tt).value == 55 % lose 25 points
%         %             lose_sample{tt} = event(tt).sample + (0: wind_samples);
%         %         end
%         %
%         %     end
%         
%         
%     end
%     
%     ratemood_sample = cell2mat(ratemood_sample); 
%     lose_sample = cell2mat(lose_sample);
%     win_sample = cell2mat(win_sample);
%     gamble_sample = cell2mat(gamble_sample);
%     certain_sample = cell2mat(certain_sample);
%     answer_sample = cell2mat(answer_sample);
%     
%     if ratemood_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
%         ratemood_sample(end,:) = [];
%     end
%     if lose_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
%         lose_sample(end,:) = [];
%     end
%     if win_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
%         win_sample(end,:) = [];
%     end
%     if gamble_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
%         gamble_sample(end,:) = [];
%     end
%     if answer_sample(end,end)> length(data_filt)  % eliminates last trial if window exceeds the length of the dataset
%         answer_sample(end,:) = [];
%     end
%     
%     if certain_sample(end,end)> length(data_filt)
%         certain_sample(end,:) = [];
%     end
%     
    
    %% Beamformer
    k = min(size(lose_sample,1),size(win_sample,1));
    
    Cc = zeros([size(C),size(lose_sample,1)]);
    for tt= 1:size(lose_sample,1)
        data_lose = data_filt(lose_sample(tt)+(0:f),:);
        Cc(:,:,tt) = cov(data_lose);
    end
    Cc = mean(Cc(:,:,randperm(size(lose_sample,1),k)),3);
%     Cc = mean(Cc,3);
    
    Ca = zeros([size(C),size(win_sample,1)]);
    for tt= 1:size(win_sample,1)
        data_win = data_filt(win_sample(tt,:),:);
        Ca(:,:,tt) = cov(data_win);
    end
%     Ca = mean(Ca,3);
    Ca = mean(Ca(:,:,randperm(size(win_sample,1),k)),3);
    
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

ttime  =0.2;
       
    zname = sprintf('PseudoZ_win_%.1fs_%d-%dHz.nii',ttime,freqband{bb}(1),freqband{bb}(2));
    
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
    
    zname = sprintf('PseudoZ_lose_%.1fs_%d-%dHz.nii',ttime,freqband{bb}(1),freqband{bb}(2));
    
%     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoZc{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        ft_write_mri(zname,sourcePostInt,'dataformat','nifti');
%     end
    
    tname =  sprintf('PseudoT_win-lose_%.1fs_%d-%dHz.nii',ttime,freqband{bb}(1),freqband{bb}(2));
        
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
    

%% Plots
close all


cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow =  PseudoZa{bb} ;

sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);

crang = [];

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcePostInt);
set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
    'position',[50 600 845 719])


cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow =  PseudoZc{bb};

sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);

crang = [];

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcePostInt);
set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
    'position',[1000 600 845 719])



cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow =  PseudoT{bb};

sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);

crang = [];

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcePostInt);
set(gcf, 'name', sprintf('Frequency %d-%dHz',freqband{bb}(1),freqband{bb}(2)),...
    'position',[1000 600 845 719])



%%

x = [81 -7 39   % right frontal 
    14 38 65     % left amygdala? 
    23 -47 64    % right amygdala
    ]   

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MRC15';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:30;
cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window
cfg.toi          = -0.5:0.05:1.5;
TFRhann7 = ft_freqanalysis(cfg, dataFIC);
To plot the result use *ft_singleplotTFR

cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.zlim         = [-3e-27 3e-27];
cfg.channel      = 'MRC15';
cfg.interactive  = 'no';
figure
ft_singleplotTFR(cfg, TFRhann7);



