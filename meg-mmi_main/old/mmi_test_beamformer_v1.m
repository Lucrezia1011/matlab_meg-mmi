clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
sub = '24138'; % with pixel '24138'; '24103';    % no pixel '24172'; '24071'
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

mri_name = [sub,'_anat_1+orig.BRIK'];
% mri_name = [sub,'_anat_orient+orig.BRIK'];

data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 1-150 Hz to adjust baseline

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
cfg.resolution      = 1;   % use a 3-D grid with a 0.5 cm resolution
cfg.sourcemodel.unit   = 'cm';
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);

%%
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

% Can either save data without frequency filtering and no baseline correction
% OR with filter and trial baseline correction
% cd(['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'])
% cfg = [];
% cfg.dataset = '24138MMI_mmi3_20190903_01.ds';
% cfg.continuous = 'yes';
% cfg.channel = 'MEG';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [1 150];
% cfg.demean ='yes';
% data = ft_preprocessing(cfg);
% 
% figure; plot(data.time{1}(1:24000), bsxfun(@plus,data.trial{1}(1:30,1:24000)', 0:1e-12:29e-12)','k')
% 
% cd(data_path)
% data_name = [sub,'MMI_mmi3_proc-f.ds']; 
% cfg = [];
% cfg.dataset = data_name;
% cfg.continuous = 'yes';
% cfg.channel = 'MEG';
% cfg.bpfilter = 'no';
% cfg.bpfreq = [0.5 180];
% cfg.demean ='no';
% data_proc = ft_preprocessing(cfg);
% 
% figure; plot(data_proc.time{1}(1:24000), bsxfun(@plus,data_proc.trial{1}(1:30,1:24000)', 0:1e-12:29e-12)','k')
% 
% n = 60;
% clf; plot(data.time{1}(12000*n+1:12000*(n+3)), data.trial{1}(29,12000*n+1:12000*(n+3)))
% hold on; plot(data_proc.time{1}(12000*n+1:12000*(n+3)), data.trial{1}(29,12000*n+1:12000*(n+3)))


ht = hilbert(data.trial{1}');
bad_segment = abs(ht)> 2e-12; % Eliminate data segments with excessive noise
bad_segment = max(bad_segment,[],2);

data.trial{1}(:,bad_segment) = [];
data.time{1}(bad_segment)= [];
data.sampleinfo = [1, length(data.time{1})];

cfg =[];
cfg.method = 'pca';
comp_pca = ft_componentanalysis(cfg, data);
score = comp_pca.trial{1}';
compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
icomp = nnz(compvar95) ;
fprintf('%d components for 95perc. of data variance\n',icomp)


% [coeff,score,latent]=pca( data.trial{1}');
% compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
% icomp = nnz(compvar95) ;
% fprintf('%d components for 95perc. of data variance\n',icomp)

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
%
% FF = fft(comp.trial{1}');
% P = abs(FF(1:length(FF)/2+1,:));
% figure; plot(linspace(0,600,length(FF)/2+1), P(:,23))

cfg           = [];
cfg.component = input('ICA component to eliminate: ');
data_clean    = ft_rejectcomponent(cfg, comp,data);


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
    ITI_match = [];
    mood_match = [];
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
            mood_match.pix_index(nn) = ii;
            mood_match.bv_index(nn) = jj;
            mood_match.sample(nn) = pix_sample(ii);
            mood_match.mood(nn) = bv.happySlider_response(jj);
        end
        
        % Rate mood
        sdiff = abs(pix_sample(ii) - (bv.blockHappySlider_started -bv_timeoffset)*f);
        [mm, jj] = min(sdiff);
        if mm < .1*f
            nn = jj;
            mood_match.pix_index(nn) = ii;
            mood_match.bv_index(nn) = jj;
            mood_match.sample(nn) = pix_sample(ii);
            mood_match.mood(nn) = bv.blockHappySlider_response(jj);
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
    
%     mm = mean(bv.outcomeAmount(outcome_match.sample>0));
%     ss = std(bv.outcomeAmount(outcome_match.sample>0));
%     
%     bigwin = (bv.outcomeAmount-mm) > ss/2;
%     bigloss = (bv.outcomeAmount-mm) < -ss/2;
%     
%     
%     for ii = 1:jj
%         win_sample = cat(2,outcome_match.sample(bigwin(1:jj) & abs(outcome_match.win)'),...
%             choice_match.sample(bigwin(1:jj) & ~outcome_match.win' ) );
%         lose_sample =cat(2,outcome_match.sample(bigloss(1:jj) & abs(outcome_match.win)'),...
%             choice_match.sample(bigloss(1:jj) & ~outcome_match.win' ) );
%     end
%     win_amount = bv.outcomeAmount(bigwin(1:jj));
%     lose_amount = bv.outcomeAmount(bigloss(1:jj));
%     
    
    pRPE = outcome_match.win == 1 ;
    nRPE = outcome_match.win == -1 ;
    pRPE_sample = outcome_match.sample(pRPE);
    nRPE_sample = outcome_match.sample(nRPE);
    
    % RPE column is the estimated RPE
    
    
    
    % the presentation of new  options (i.e. answer match) should induce a
    % signal representing the probability of winning, based on current mood
    % and values of options, influencing behaviour (gamble/certain choice).
    
    % How is mood represented?? Change of electrical potential of a set of
    % neurons? Could it change the baseline power in specific bands.
    
end


%%

freqband =cell(1,4);
freqband{1} = [1,4];
freqband{2} = [4,8];
freqband{3} = [8,13];
freqband{4} = [13,30];

for bb=1:4
    filt_order = [200]; % default
    data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, freqband{bb},filt_order,'firls')';
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
        nfft = 2^13;
        F = fft(data_filt, nfft);
        P = abs(F(1:nfft/2+1, :));
        freqs = linspace(0,data.fsample/2, nfft/2+1);
        figure(1); clf; plot(freqs, mean(P,2))
    
    
    %% Beamformer
    
    
    % Cue-P3 ERP:
    
    wind = 0*f:0.6*f;
    % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    
    win_logic = bv.RPE > 10 & bv.outcomeAmount > 10;
    win_logic = bv.RPE > 10;
    win_logic = win_logic(1:length(outcome_match.sample));
    win_samples = outcome_match.sample(win_logic);
    
    
    lose_logic = bv.RPE < -10 & bv.outcomeAmount < -10;
    lose_logic = bv.RPE < -10;
    lose_logic = lose_logic(1:length(outcome_match.sample));
    lose_samples = outcome_match.sample(lose_logic);
    
    
    
    Cap = zeros([size(C),size(win_samples,2)]);
    for tt= 1:size(win_samples,2)        
        data_win = data_filt(win_samples(tt)+wind,:);
        Cap(:,:,ii) = cov(data_win);
    end   
    Cap = mean(Cap,3);
    
    
    Can = zeros([size(C),size(lose_samples,2)]);
    for tt= 1:size(lose_samples,2)        
        data_lose = data_filt(lose_samples(tt)+wind,:);
        Can(:,:,ii) = cov(data_lose);
    end
    Can = mean(Can,3);
    

%     ii = 0;
%     Cap = zeros([size(C),nnz(answer_match.sample)]);
%     for tt= 1:size(answer_match.sample,2)
%         if answer_match.sample(tt) > 0
%             ii = ii+1;
%             data_win = data_filt(answer_match.sample(tt)+wind,:);
%             Cap(:,:,ii) = cov(data_win);
%         end
%     end
%     
%     Cap = mean(Cap,3);
    
    % Cue-P3 ERP:
    
    wind = 0.7*f:1.3*f;
    ii = 0;
    Cc = zeros([size(C),nnz(answer_match.sample)]);
    for tt= 1:size(ITI_match.sample,2)
        if ITI_match.sample(tt) > 0
            ii = ii +1;
            data_win = data_filt(ITI_match.sample(tt)+wind,:);
            Cc(:,:,ii) = cov(data_win);
        end
    end
    
    if nnz(answer_match.sample)> nnz(ITI_match.sample)
        for tt = ii+1:nnz(answer_match.sample)
            data_win = data_filt(ITI_match.sample(end)+wind-0.3*f,:);
            Cc(:,:,tt) = cov(data_win);
        end
    elseif nnz(answer_match.sample)< nnz(ITI_match.sample)
        Cc = Cc(:,:,1:nnz(answer_match.sample));
    end
    Cc= mean(Cc,3);
    
    
%     Cna = eye(size(Ca))*min(svd(Ca));
%     Cnc = eye(size(Cc))*min(svd(Cc));
    %% Can make timecourse of covariance
    L = grid.leadfield;
    W = cell(size(L));
    
    dipout = cell(size(L));
    
    Za(1:size(L,2)) = {0};
    Zc(1:size(L,2)) = {0};
    Tstat(1:size(L,2)) = {0};
    Tstatn(1:size(L,2)) = {0};
    
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
            
            
            %             na = w'*Cna*w;
            %             nc = w'*Cnc*w;
            
            na = w'*noise*w;
            nc = w'*noise*w;
            
            Qap = w'*Cap*w;
            Qan = w'*Can*w;
            Qc = w'*Cc*w;
            
            
            W{ii} = w;
%             Tstat{ii} = (Qa - Qc) ./ (na + nc);
            Tstatn{ii} = (Qap - Qan) ./ (Qa + Qc); % normalised version
            
            Za{ii} = Qap / na;
            Zc{ii} = Qan / na;
            
            %
            %             wnorm = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo );
            %             Qa = wnorm'*Ca*wnorm;
            %             Qc = wnorm'*Cc*wnorm;
            
            %             Za{ii} = Qa;
            %             Zc{ii} = Qc;
        end
    end
%     PseudoT{bb} = cell2mat(Tstat);
    PseudoTn{bb} = cell2mat(Tstatn);
    PseudoZa{bb} = cell2mat(Za);
    PseudoZc{bb} = cell2mat(Zc);
    
     cfg = [];
        cfg.parameter = 'pow';
        sourceTstat = struct;
        sourceTstat.dim = grid.dim;
        sourceTstat.inside = grid.inside;
        sourceTstat.pos = grid.pos;
        sourceTstat.method = 'average';
        
        %     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoTn{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname =  sprintf('RPEdiff_0.6s_PseudoT_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
        zname_nii = [zname,'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
    
    
        sourceTstat.avg.pow =  PseudoZa{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname =  sprintf('RPEpos_0.6s_PseudoZ_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
        zname_nii = [zname,'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
        
        
        sourceTstat.avg.pow =  PseudoZc{bb};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname =  sprintf('RPEneg_0.6s_PseudoZ_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
        zname_nii = [zname,'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
        
     
    %% Time 0 to 1050ms (300ms window, 150ms step) following outcome with Positive RPE
    % keyboard
    tstep = 0.15;
    twind = 0.3;
    
    wind = 1*f:(1+twind)*f;
    Ccall = zeros([size(C),nnz(ITI_match.sample)]);
    ii = 0;
    for tt= 1:length(ITI_match.sample)
        if ITI_match.sample(tt) > 0
            ii = ii+1;
            data_win = data_filt(ITI_match.sample(tt)+wind,:);
            Ccall(:,:,ii) = cov(data_win);
        end
    end
    Ccall= mean(Ccall,3);
    
    
    % Cue-P3 ERP:
    % Can make 4D (sliding time window) Tstat of:
    % Answer to choice (separate certain and gamble), can do 0 to 600ms and
    % -600ms from button press
    
    % Choice to outcome (separate gamble and certain)
    
    % Outcome (gamble), separate positive and negative RPE.
    
    n_wind = 6;
    
    time = zeros(2,n_wind);
    Ca4D = zeros([size(C),n_wind]);
    % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    
   
    
    na = zeros(1,n_wind);
    for tt4 = 1:n_wind
        wind = ((tt4-1)*tstep*f) + (0:twind*f);
        time(:,tt4) = ((tt4-1)*tstep) + [0,twind];
        
        Ca = zeros([size(C),size(pRPE_sample,2)]);
        for tt= 1:size(pRPE_sample,2)
            
            data_win = data_filt(pRPE_sample(tt)+wind,:);
            Ca(:,:,tt) = cov(data_win);
        end
        
        Ca4D(:,:,tt4) = mean(Ca,3);
        na(tt4) = min(svd(Ca4D(:,:,tt4)));
        
    end
    
    na = mean(na).*eye(size(C));
    
    
    % Cue-P3 ERP:
    
    wind = 1*f:(1+twind)*f;
    Cc = zeros([size(C),size(pRPE_sample,2)]);
    
    ITI_sample = ITI_match.sample(pRPE);
    for tt= 1:size(pRPE_sample,2)
        data_win = data_filt(ITI_sample(tt)+wind,:);
        Cc(:,:,tt) = cov(data_win);
        
    end
    Cc= mean(Cc,3);
    
    nc = min(svd(Cc)).*eye(size(C));
    
    
    PseudoT4D = cell(1,n_wind);
    for tt = 1:n_wind
        Tstat(1:size(L,2)) = {0};
        for ii = 1:length(W)
            w = W{ii};
            if ~isempty(w)
                Qa = w'*Ca4D(:,:,tt)*w;
                %             Qc = w'*Cc*w;
                Qc = w'*Ccall*w;
%                 Tstat{ii} = (Qa - Qc) ./ (2*w'*noise*w);
                Tstat{ii} = (Qa - Qc) ./ (2*Qc);
            end
        end
        PseudoT4D{tt} = cell2mat(Tstat);
    end
    
    
    ttime  =length(wind)/f;
    
    cd(data_path)
    
    if ~exist('result','dir')
        mkdir('result')
    end
    
    cd result
    
    zname =  sprintf('PositiveRPE_PseudoT_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
    zname_nii_cat = [];
    
    for tt = 1:n_wind
        %     zname =  sprintf('PositiveRPE_PseudoT_%.1fs_%d-%dHz_%d',ttime,freqband{bb}(1),freqband{bb}(2),tt);
        
        cfg = [];
        cfg.parameter = 'pow';
        sourceTstat = struct;
        sourceTstat.dim = grid.dim;
        sourceTstat.inside = grid.inside;
        sourceTstat.pos = grid.pos;
        sourceTstat.method = 'average';
        
        %     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoT4D{tt};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname_nii = [zname,'_',num2str(tt),'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
        
        zname_nii_cat = cat(2,zname_nii_cat,' ',zname_nii);
        
        %     unix(['3dcopy ',zname_nii,' ',[zname,'+orig',num2str(tt),'.BRIK']])
    end
    
    unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_tcat ',zname_nii_cat]);
%     unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_fixnoise_tcat ',zname_nii_cat]);
    unix(['rm ',zname_nii_cat])
    
    
    %% Time 0 to 1050ms (300ms window, 150ms step) following outcome with Negative RPE
    % keyboard
    
    
    % Cue-P3 ERP:
    % Can make 4D (sliding time window) Tstat of:
    % Answer to choice (separate certain and gamble), can do 0 to 600ms and
    % -600ms from button press
    
    % Choice to outcome (separate gamble and certain)
    
    % Outcome (gamble), separate positive and negative RPE.
    
    n_wind = 6;
    
    time = zeros(2,n_wind);
    Ca4D = zeros([size(C),n_wind]);
    % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    
    tstep = 0.15;
    twind = 0.3;
    
    na = zeros(1,n_wind);
    
    for tt4 = 1:n_wind
        wind = ((tt4-1)*tstep*f) + (0:twind*f);
        time(:,tt4) = ((tt4-1)*tstep) + [0,twind];
        
        Ca = zeros([size(C),size(nRPE_sample,2)]);
        for tt= 1:size(nRPE_sample,2)
            
            data_win = data_filt(nRPE_sample(tt)+wind,:);
            Ca(:,:,tt) = cov(data_win);
        end
        
        Ca4D(:,:,tt4) = mean(Ca,3);
        na(tt4) = min(svd(Ca4D(:,:,tt4)));
    end
    
    na = mean(na).*eye(size(C));
    
    % Cue-P3 ERP:
    
    wind = 1*f:(1+twind)*f;
    Cc = zeros([size(C),size(nRPE_sample,2)]);
    
    ITI_sample = ITI_match.sample(pRPE);
    for tt= 1:size(nRPE_sample,2)
        data_win = data_filt(ITI_sample(tt)+wind,:);
        Cc(:,:,tt) = cov(data_win);
        
    end
    Cc= mean(Cc,3);
    nc = min(svd(Cc)).*eye(size(C));
    
    
    
    PseudoT4D = cell(1,n_wind);
    for tt = 1:n_wind
        Tstat(1:size(L,2)) = {0};
        for ii = 1:length(W)
            w = W{ii};
            if ~isempty(w)
                Qa = w'*Ca4D(:,:,tt)*w;
                %             Qc = w'*Cc*w;
                Qc = w'*Ccall*w;
%                 Tstat{ii} = (Qa - Qc) ./ (2*w'*noise*w);
                Tstat{ii} = (Qa - Qc) ./ (2*Qc);
            end
        end
        PseudoT4D{tt} = cell2mat(Tstat);
    end
    
    
    ttime  =length(wind)/f;
    
    cd(data_path)
    
    if ~exist('result','dir')
        mkdir('result')
    end
    
    cd result
    
    zname =  sprintf('NegativeRPE_PseudoT_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
    zname_nii_cat = [];
    
    for tt = 1:n_wind
        %     zname =  sprintf('PositiveRPE_PseudoT_%.1fs_%d-%dHz_%d',ttime,freqband{bb}(1),freqband{bb}(2),tt);
        
        cfg = [];
        cfg.parameter = 'pow';
        sourceTstat = struct;
        sourceTstat.dim = grid.dim;
        sourceTstat.inside = grid.inside;
        sourceTstat.pos = grid.pos;
        sourceTstat.method = 'average';
        
        %     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoT4D{tt};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname_nii = [zname,'_',num2str(tt),'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
        
        zname_nii_cat = cat(2,zname_nii_cat,' ',zname_nii);
        
        %     unix(['3dcopy ',zname_nii,' ',[zname,'+orig',num2str(tt),'.BRIK']])
    end
    
    unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_tcat ',zname_nii_cat]);
%     unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_fixnoise_tcat ',zname_nii_cat]);
    unix(['rm ',zname_nii_cat])

    
    %% Positive vs Negative RPE
    % keyboard    
    
    % Cue-P3 ERP:
    % Can make 4D (sliding time window) Tstat of:
    % Answer to choice (separate certain and gamble), can do 0 to 600ms and
    % -600ms from button press
    
    % Choice to outcome (separate gamble and certain)
    
    % Outcome (gamble), separate positive and negative RPE.
    
    n_wind = 6;
    
    time = zeros(2,n_wind);
    Ca4D = zeros([size(C),n_wind]);
    Cc4D = zeros([size(C),n_wind]);
    % cuep3_sample = zeros(nnz(choice_match.sample), length(wind));
    
    tstep = 0.15;
    twind = 0.3;
    
    na = zeros(1,n_wind);
    
    for tt4 = 1:n_wind
        wind = ((tt4-1)*tstep*f) + (0:twind*f);
        time(:,tt4) = ((tt4-1)*tstep) + [0,twind];
        
        Ca = zeros([size(C),size(pRPE_sample,2)]);
        for tt= 1:size(pRPE_sample,2)
            
            data_win = data_filt(pRPE_sample(tt)+wind,:);
            Ca(:,:,tt) = cov(data_win);
        end
        
        Ca4D(:,:,tt4) = mean(Ca,3);
        na(tt4) = min(svd(Ca4D(:,:,tt4)));
    end
    
    na = mean(na).*eye(size(C));
    
   nc = zeros(1,n_wind);
    for tt4 = 1:n_wind
        wind = ((tt4-1)*tstep*f) + (0:twind*f);
        time(:,tt4) = ((tt4-1)*tstep) + [0,twind];
        
        Cc = zeros([size(C),size(nRPE_sample,2)]);
        for tt= 1:size(nRPE_sample,2)
            
            data_win = data_filt(nRPE_sample(tt)+wind,:);
            Cc(:,:,tt) = cov(data_win);
        end
        
        Cc4D(:,:,tt4) = mean(Cc,3);
        nc(tt4) = min(svd(Cc4D(:,:,tt4)));
    end
    
    nc = mean(nc).*eye(size(C));
      
    
    PseudoT4D = cell(1,n_wind);
    PseudoZp = cell(1,n_wind);
    PseudoZn = cell(1,n_wind);
    for tt = 1:n_wind
        Tstat(1:size(L,2)) = {0};
        Za(1:size(L,2)) = {0};
        Zc(1:size(L,2)) = {0};
        for ii = 1:length(W)
            w = W{ii};
            if ~isempty(w)
                Qa = w'*Ca4D(:,:,tt)*w;
                %             Qc = w'*Cc*w;
                Qc = w'*Cc4D(:,:,tt)*w;
%                 Tstat{ii} = (Qa - Qc) ./ (2*w'*noise*w);
                Tstat{ii} = (Qa - Qc) ./ (Qa + Qc);
                Za{ii} = Qa / (w'*na*w);
                Zc{ii} = Qc / (w'*nc*w);
            end
        end
        PseudoT4D{tt} = cell2mat(Tstat);
        PseudoZp{tt} = cell2mat(Za);
        PseudoZn{tt} = cell2mat(Zc);
    end
    
    
    ttime  =length(wind)/f;
    
    cd(data_path)
    
    if ~exist('result','dir')
        mkdir('result')
    end
    
    cd result
    
    zname =  sprintf('RPEdiff_PseudoT_%d-%dHz',freqband{bb}(1),freqband{bb}(2));
    zname_nii_cat = [];
    
    for tt = 1:n_wind
        %     zname =  sprintf('PositiveRPE_PseudoT_%.1fs_%d-%dHz_%d',ttime,freqband{bb}(1),freqband{bb}(2),tt);
        
        cfg = [];
        cfg.parameter = 'pow';
        sourceTstat = struct;
        sourceTstat.dim = grid.dim;
        sourceTstat.inside = grid.inside;
        sourceTstat.pos = grid.pos;
        sourceTstat.method = 'average';
        
        %     if ~exist(zname,'file')
        
        sourceTstat.avg.pow =  PseudoT4D{tt};
        sourcePostInt  = ft_sourceinterpolate(cfg, sourceTstat , mri);
        sourcePostInt.anatomy = sourcePostInt.pow;
        
        zname_nii = [zname,'_',num2str(tt),'.nii'];
        ft_write_mri(zname_nii,sourcePostInt,'dataformat','nifti');
        
        zname_nii_cat = cat(2,zname_nii_cat,' ',zname_nii);
        
        %     unix(['3dcopy ',zname_nii,' ',[zname,'+orig',num2str(tt),'.BRIK']])
    end
    
    unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_tcat ',zname_nii_cat]);
%     unix(['3dTcat -tpattern seqplus -tr ',num2str(tstep),' -prefix ',zname,'_fixnoise_tcat ',zname_nii_cat]);
    unix(['rm ',zname_nii_cat])
    
    %%
end

%% Compare mood

plot(mood_match.sample(mood_match.sample>0)/f, mood_match.mood(mood_match.sample>0))

x = mood_match.sample(mood_match.sample>0)';
x(:,2) = 1;
v = mood_match.mood(mood_match.sample>0)';
xi = x(1):length(data_clean.time{1});
xi(2,:) = 1;
vq = griddatan(x,v,xi');

x = mood_match.sample(mood_match.sample>0)';
v = mood_match.mood(mood_match.sample>0)';
F = griddedInterpolant(x,v,'pchip');
xi = x(1):length(data_clean.time{1});
vi = F(xi); % Mood timecourse

hold on
plot(xi,vi)