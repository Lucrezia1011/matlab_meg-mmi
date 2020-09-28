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

% Check timecourses of run1 and run2 have comparable amplitude!!
for sn = [16] % redo 15, should keep same number of sensors % all subjects with continuos recordings
clearvars -except sn subn
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

%%
for runs = 1:length(data_names)
    h = ft_read_header(data_names{runs});
    nsamples(runs) = h.nTrials;
end
[~,runorder]=sort(nsamples); % start from shorter recording
nsamples = 1200*20*60;
noiseC = [];
[t,mood] = plot_mood(sub,false,false);
close

Pruns = [];
for runs = runorder
data_name = data_names{runs};

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
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

delete([processing_folder,'/lead_fields.mat'])
%% Clean data
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
% cfg.bpfreq = [110 140];
% cfg.bpfilter = 'yes';
data = ft_preprocessing(cfg);

filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [110 140],filt_order,'but');

zchans = zscore(data_filt');

for ii = 1:size(zchans,2)
    zchans(:,ii) = smooth(abs(zchans(:,ii)),10);
end

tempchans = false(length(data.label),1);
for ii = 1:length(data.label)
    tempchans(ii) = strncmp(data.label{ii},'MRT',3) || strncmp(data.label{ii},'MLT',3);
end
% musclez = mean(zchans(:,tempchans)>3,2);
musclez = mean(zchans>3,2);

% figure; plot(data.time{1},musclez) %zscore >3 in temporal channels
ind = find(musclez>0.3); %30% of all channels channels have artifacts
muscleind = cell(1,length(ind));
for ii = 1:length(ind)
    muscleind{ii} = (ind(ii)-data.fsample) : (ind(ii)+data.fsample); % pad by 1s 
end
muscleind = unique(cell2mat(muscleind));
    
% hold on
% plot(muscleind/data.fsample,ones(1,length(muscleind))*0.4,'rx')



% 
% cfg = [];
% cfg.dataset = data_name;
% cfg.continuous = 'yes';
% % cfg.trl = cat(2, data.sampleinfo, -ones(size(data.sampleinfo,1),1)*240);
% % cfg.trl(:,1) = cfg.trl(:,1)+240;
% % cfg.trl(:,2) = cfg.trl(:,2)-240;
% cfg.trl = [241 data.sampleinfo(end)-240 -240]; 
% [cfg, artifact] = ft_artifact_muscle(cfg);
% 
% trial = 75;
% figure(1); clf 
% plot(data.sampleinfo(trial,1):data.sampleinfo(trial,2), (data.trial{trial}(98:131,:)'),'k')
% hold on
% plot(data.sampleinfo(trial,1):data.sampleinfo(trial,2), ( data.trial{trial}(227:260,:)')- 4e-13,'k')
% ind = true ;
% ii = 1;
% while ind
%     if artifact(ii,2)<=data.sampleinfo(trial,2) && artifact(ii,1)>=data.sampleinfo(trial,1) 
%     plot(artifact(ii,1):artifact(ii,2),-ones(1,(artifact(ii,2)-artifact(ii,1))+1)*2e-13, 'r','linewidth',2)    
%     elseif artifact(ii,2)>data.sampleinfo(trial,2)
%         ind = false;
%     end
%     ii=ii+1;
% end
% ylim([-8 4]*1e-13)


f = data.fsample;

if exist([processing_folder,'/ICA_artifacts.mat'],'file')
    load([processing_folder,'/ICA_artifacts.mat']);
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data_clean    = ft_rejectcomponent(cfg, comps,data);
end

data = data_clean;
clear data_clean

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
%%
bvf = fieldnames(bv_match);

cellindc = [];
for jj = 2:(length(bvf)-2) %exclude time, ITI and buttonpress
    ind = bv_match.(bvf{jj}).sample; 
    ind(ind==0) = [];
    cellind = cell(1,length(ind));
    for ii = 1:length(ind)
        cellind{ii} = (ind(ii)-data.fsample*0.1) : (ind(ii)+data.fsample*0.5); % pad by 1s 
    end
    cellindc = cat(2,cellindc, cellind);
end
cellindc = unique(cell2mat(cellindc));


ind = buttonpress;
cellind = cell(1,length(buttonpress));
for ii = 1:length(ind)
    cellind{ii} = (ind(ii)-data.fsample*0.1) : (ind(ii)+data.fsample*0.5); % pad by 1s 
end
buttonind = unique(cell2mat(cellind));
 
inddel = unique([cellindc,buttonind,muscleind]);
inddel(inddel<1) = [];
inddel(inddel>data.sampleinfo(2)) = [];
indkeep = 1:data.sampleinfo(end);
indkeep = setdiff(indkeep,inddel);
% da
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, [1 4],200,'firls');
% 
% 
% 
% trials = blockmood_match.sample;
% trials(trials==0) = [];
% trial1 = data.trial{1}(:,intersect(trials(1):trials(2),indkeep));
% trial2 = data.trial{1}(:,intersect(trials(2):trials(3),indkeep));
% trial3 = data.trial{1}(:,intersect(trials(3):data.sampleinfo(end),indkeep));
% 
% nn = 2^16;
% F1 = fft(trial1',nn);
% figure; plot(linspace(0,f,nn),abs(F1)); xlim([0 200])
% 
% F2 = fft(trial2',nn);
% figure; plot(linspace(0,f,nn),abs(F2)); xlim([0 200])
% 
% F3 = fft(trial3',nn);
% figure; plot(linspace(0,f,nn),abs(F3)); xlim([0 200])

%% TFS
% data_cut = data;
% data_cut.time{1}(inddel) = [];
% data_cut.trial{1}(:,inddel) = [];
% data_cut.sampleinfo(2) = length(indkeep);
% 
% cfg = [];
% cfg.length  = 1;
% cfg.overlap = 0;
% data_segmented = ft_redefinetrial(cfg, data_cut);
% 
% cfg = [];
% cfg.method     = 'mtmfft';
% cfg.taper      = 'hanning';
% cfg.foilim     = [1 40];
% cfg.keeptrials = 'yes';
% freq_segmented = ft_freqanalysis(cfg, data_segmented);
% begsample = data_segmented.sampleinfo(:,1);
% endsample = data_segmented.sampleinfo(:,2);
% time = ((begsample+endsample)/2) / data_segmented.fsample;
% 
% figure; 
% F = squeeze(freq_segmented.powspctrm(:,100,:));
% F = zscore(F);
% pcolor(time,freq_segmented.freq, F')
% shading interp; caxis([-3 3])
% 
% freq_continuous           = freq_segmented;
% freq_continuous.powspctrm = permute(freq_segmented.powspctrm, [2, 3, 1]);
% freq_continuous.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
% freq_continuous.time      = time;  
% cfg = [];
% cfg.layout = 'CTF275_helmet.mat';
% cfg.zlim  = [0 1e-26];
% ft_multiplotTFR(cfg, freq_continuous);
%% Beamfomer

data_cut = data;
data_cut.time{1}(inddel) = [];
data_cut.trial{1}(:,inddel) = [];
data_cut.sampleinfo(2) = length(indkeep);

data_filt = ft_preproc_bandpassfilter(data_cut.trial{1}, data.fsample, [1 150],[],'but');
data_cut.trial{1} = data_filt;

cfg = [];
cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
data_cut = ft_resampledata(cfg, data_cut);

% Need to adjust time to reflect cut points
data_cut.time{1} = data.time{1};
data_cut.time{1}(inddel) = [];
data_cut.time{1} = downsample(data_cut.time{1},4);


icacomps = length(data.cfg.component);
if length(data_cut.trial{1}) < nsamples
    nsamples = length(data_cut.trial{1});
end
C = cov(data_cut.trial{1}(:,1:nsamples)');

if isempty(noiseC)
    E = svd(C);
    nchans = length(data.label);
    noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components
end
Cr = C + 4*noiseC; % need to normalise because of ICA
% Cr = C + 0.01*eye(nchans)*E(1);
L = grid.leadfield;
% 
% % VEp(1:size(L,2)) = {0};
% % VEn(1:size(L,2)) = {0};
% 
% % VE = cell(1,size(L,2));
% W(1:size(L,2)) = {0};
% 
% %%
% 
% P = [];
% nv = 200;
% VE = zeros(nv,data_cut.sampleinfo(2));
% 
% data_ve = data_cut;
% data_ve.label = data_ve.label(1:nv);
% data_ve.trial = [];
% 
% kk = 0;
% 
% for ii = 1:length(L)
%     lf = L{ii}; % Unit 1Am
% 
%     if ~isempty(lf)
%         
%             kk = kk+1;
%             % %  G O'Neill method, equivalent to ft
%             [v,d] = svd(lf'/Cr*lf);
%             d = diag(d);
%             jj = 2;
%             %         if d(3) < 1
%             %             jj = 2; % The minumum singular value is degenerate
%             %         else
%             %             jj =3;
%             %         end
%             lfo = lf*v(:,jj); % Lead field with selected orientation
%             w = Cr\lfo / (lfo'/Cr*lfo) ;       
%             %         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
%             %         Better Hall's or normalized weights
%             wnorm = w/sqrt(w'*noiseC*w);
%             VE(kk,:) = wnorm'*data_filt;     
%             
%             if kk == nv
%                 
%                 clc
%                 data_ve.trial{1} = VE;
%                 cfg = [];
%                 cfg.length  = 10;
%                 cfg.overlap = 0.8;
%                 data_segmented = ft_redefinetrial(cfg, data_ve);
% 
%                 cfg = [];
%                 cfg.method     = 'mtmfft';
%                 cfg.taper      = 'hanning';
%                 cfg.pad        ='nextpow2';
%                 cfg.foi        = [1:58];
%                 cfg.keeptrials = 'yes';
%                 freq_segmented = ft_freqanalysis(cfg, data_segmented);
%                 P = cat(2,P,freq_segmented.powspctrm);
%                 kk=0;
%                 fprintf('Done %.0f perc.\n', ii/length(L)*100)
%                 VE = zeros(nv,data_cut.sampleinfo(2));
% 
%             end
% %              
%     end   
%     
% end
% 
% clc
% data_ve.trial{1} = VE(1:kk,:);
% data_ve.label = data_ve.label(1:kk);
% cfg = [];
% cfg.length  = 10;
% cfg.overlap = 0.8;
% data_segmented = ft_redefinetrial(cfg, data_ve);
% 
% cfg = [];
% cfg.method     = 'mtmfft';
% cfg.taper      = 'hanning';
% cfg.pad        ='nextpow2';
% cfg.foi        = [1:58];
% cfg.keeptrials = 'yes';
% freq_segmented = ft_freqanalysis(cfg, data_segmented);
% P = cat(2,P,freq_segmented.powspctrm);
% kk=0;
% fprintf('Done %.0f perc.\n', ii/length(L)*100)

%% Power spectra for each mood rating

nfft = 2048;
P = cell(size(L));

k = 0;
parfor ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)
            k = k+1;
            % %  G O'Neill method, equivalent to ft
            [v,d] = svd(lf'/Cr*lf);
            d = diag(d);
            jj = 2;
           
            lfo = lf*v(:,jj); % Lead field with selected orientation
            w = Cr\lfo / (lfo'/Cr*lfo) ;       
            %         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
            %         Better Hall's or normalized weights
            wnorm = w/sqrt(w'*noiseC*w);
            VE = wnorm'*data_cut.trial{1};   
            
            Pxx = zeros(length(t{runs}),nfft/2+1);
            for tt = 1:length(t{runs})
                tm = t{runs}(tt); % mood time
                ind = (data_cut.time{1} > (tm-15)) & (data_cut.time{1} < (tm+15)); % 30s window around time of mood rating
                data_mood = VE(:,ind);
                Fxx  = fft(data_mood',nfft);
                Pxx(tt,:) = abs(Fxx(1:(nfft/2+1),:));
%                   [Pxx(tt,:),F] = pwelch(data_mood',[],[],nttf,data_cut.fsample);
            end
            P{ii} = Pxx;
           
    end
    if mod(ii,300) == 0
        clc
        fprintf('Done %.0f perc.\n', ii/length(L)*100)
    end
end
F = linspace(0,data_cut.fsample/2,nfft/2+1);
Pruns = cat(1,Pruns,P);
clear P
end

%% Plots
% k = 2016;
% figure(1); clf
% subplot(231); pcolor(t{1},F(1:(nttf/8)),P{k}(:,1:(nttf/8))')
% shading interp
% xlabel('Time (s)'); ylabel('Frequency (Hz)')
% caxis([0 600]); title('TFS from fft')
% 
% subplot(232); pcolor(t{1},F((nttf/8):end),P{k}(:,(nttf/8):end)')
% shading interp
% xlabel('Time (s)'); ylabel('Frequency (Hz)')
% caxis([0 200]); title('TFS from fft')
% 
% 
% subplot(234); pcolor(t{1},F(1:(nttf/8)),Pw{k}(:,1:(nttf/8))')
% shading interp
% xlabel('Time (s)'); ylabel('Frequency (Hz)')
% caxis([0 1]); title('TFS from pwelch')
% 
% subplot(235); pcolor(t{1},F((nttf/8):end),Pw{k}(:,(nttf/8):end)')
% shading interp
% xlabel('Time (s)'); ylabel('Frequency (Hz)')
% caxis([0 0.06]); title('TFS from pwelch')
% 
% subplot(233)
% hold all
% bar(1, corr(mean(P{k}(:,F>=1 & F<=4),2),mood{1}) )
% bar(2, corr(mean(P{k}(:,F>=4 & F<=8),2),mood{1}))
% bar(3, corr(mean(P{k}(:,F>=8 & F<=13),2),mood{1}) )
% bar(4, corr(mean(P{k}(:,F>=13 & F<=30),2),mood{1}) )
% bar(5, corr(mean(P{k}(:,F>=30 & F<=55),2),mood{1}))
% bar(6, corr(mean(P{k}(:,F>=65 & F<=115),2),mood{1}) )
% ylabel('Pearson correlation (r)')
% set(gca,'XTick',1:6,'XTickLabel',{'\delta';'\theta';'\alpha';'\beta';'\gamma_{low}';'\gamma_{high}'})
% ylim([-1 1])
% title('Power spectrum correlation with mood')
% 
% subplot(236)
% hold all
% bar(1, corr(mean(Pw{k}(:,F>=1 & F<=4),2),mood{1}) )
% bar(2, corr(mean(Pw{k}(:,F>=4 & F<=8),2),mood{1}))
% bar(3, corr(mean(Pw{k}(:,F>=8 & F<=13),2),mood{1}) )
% bar(4, corr(mean(Pw{k}(:,F>=13 & F<=30),2),mood{1}) )
% bar(5, corr(mean(Pw{k}(:,F>=30 & F<=55),2),mood{1}))
% bar(6, corr(mean(Pw{k}(:,F>=65 & F<=115),2),mood{1}) )
% ylabel('Pearson correlation (r)')
% set(gca,'XTick',1:6,'XTickLabel',{'\delta';'\theta';'\alpha';'\beta';'\gamma_{low}';'\gamma_{high}'})
% ylim([-1 1])
% title('Power spectrum correlation with mood')

%% 
delete Mood_corr_*

k = 0;
R = zeros(6,nnz(grid.inside));
mood = cell2mat(mood');

for ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)
        k = k+1;
        P = cell2mat(Pruns(:,ii));
        
        R(1,k) = corr(mean(P(:,F>=1 & F<=4),2),mood);
        R(2,k) = corr(mean(P(:,F>=4 & F<=8),2),mood);
        R(3,k) = corr(mean(P(:,F>=8 & F<=13),2),mood);
        R(4,k) = corr(mean(P(:,F>=13 & F<=30),2),mood);
        R(5,k) = corr(mean(P(:,F>=30 & F<=55),2),mood);
        R(6,k) = corr(mean(P(:,F>=65 & F<=115),2),mood);
    end
    
end

freqnames = {'delta';'theta';'alpha';'beta';'gamma_low';'gamma_high'};
for k = 1:6
    Ranat = zeros(1,length(L));
    Ranat(grid.inside) = R(k,:);

    sourceant =[];
    sourceant.pow = Ranat;
    sourceant.dim = grid.dim;
    sourceant.inside = grid.inside;
    sourceant.pos = grid.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
    sourceant_Int.anatomy = sourceant_Int.pow;
    writebrik(['Mood_corr_',freqnames{k}],sourceant_Int)
end
% crang = [0.2 0.5];
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourceant_Int);
% 
% crang = [-0.5 -0.2];
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourceant_Int);

end
%%

% 
% begsample = data_segmented.sampleinfo(:,1);
% endsample = data_segmented.sampleinfo(:,2);
% % time = ((begsample+endsample)/2) / data_segmented.fsample;
% 
% time = round((begsample+endsample)/2);
% time = indkeep(time)/f;
% figure; 
% x = mean(P(:,:,1:4),3)';
% pcolor(time,1:size(P,2),x)
% shading interp; 
% caxis([0 8]); title('delta')
% 
% figure; 
% x = mean(P(:,:,4:8),3)';
% pcolor(time,1:size(P,2),x)
% shading interp; 
% caxis([0 8]); title('theta')
% 
% figure; 
% x = mean(P(:,:,8:13),3)';
% pcolor(time,1:size(P,2),x)
% shading interp; 
% caxis([0 5]); title('alpha')
% 
% figure; 
% x = mean(P(:,:,13:30),3)';
% pcolor(time,1:size(P,2),x)
% shading interp; 
% caxis([0 0.5]); title('beta')
% [t,mood] = plot_mood(sub,false,false);
% close 
% F = griddedInterpolant(t{1},mood{1},'pchip');
% moodint = F(time);
% 
% X = [mean(P(:,:,1:4),3), mean(P(:,:,4:8),3), mean(P(:,:,8:13),3), mean(P(:,:,13:30),3)];
% 
% 
% X =  mean(P(:,:,1:4),3);
% 
% pvalue = cell(1,size(X,2));
% beta = cell(1,size(X,2)); 
% parfor ii= 1:size(X,2)
%     Xs = smooth(X(:,ii));
%     mdl = fitglm(moodint',Xs,'linear','Distribution','normal'); 
%     beta{ii} = mdl.Coefficients.Estimate(2);
%     pvalue{ii} = mdl.Coefficients.pValue(2);
% end
% 
% beta = cell2mat(beta);
% pvalue = cell2mat(pvalue);
% 
% [psort,ind]=sort(pvalue); 
% indsig = ind(psort<= (0.05./(length(pvalue):-1:1)));
% psort = zeros(size(pvalue));
% psort(indsig) = 1;
% 
% pvaluesig = NaN(size(grid.inside));
% pvaluesig(grid.inside) = psort;
% % ind = find(pvaluesig>(0.05/nnz(grid.inside)));
% ind = find(pvaluesig~=1);
% betasig = NaN(size(grid.inside));
% betasig(grid.inside) = beta;
% betasig(ind) = 0;
% 
% sourceant =[];
% sourceant.pow = betasig;
% sourceant.dim = grid.dim;
% sourceant.inside = grid.inside;
% sourceant.pos = grid.pos;
% cfg = [];
% cfg.parameter = 'pow';
% sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
% 
% 
% crang = [];
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourceant_Int);
% 
% 
% 
% 
