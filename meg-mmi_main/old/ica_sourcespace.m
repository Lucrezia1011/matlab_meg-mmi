clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

%% Co-register MRI from fiducial positions
sub = '24138'; % with pixel '24138'artefacts; '24103';    % no pixel '24172'; '24071'
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

%% Artifact removal
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 120];
cfg.demean = 'yes';
data = ft_preprocessing(cfg);

cfg =[];
cfg.method = 'pca';
comp_pca = ft_componentanalysis(cfg, data);
score = comp_pca.trial{1}';
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
cfg.channel  = [1:5]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, comp)

clc
cfg           = [];
cfg.component = input('ICA component to eliminate: ');
data_clean    = ft_rejectcomponent(cfg, comp,data);

keyboard
%% Beamforming
% data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [1 55], 200, 'firls', 'twopass')';
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [0.5 55], [], 'but', 'twopass')';

mu = 10;
C = cov(data_filt);
noise = min(svd(C)).*eye(size(C));
% Cr = C + mu*noise;
Cr = C + 0.01*max(svd(C)).*eye(size(C));

L = grid.leadfield;
W = cell(size(L));
VE = cell(1,size(L,2));

parfor ii = 1:length(L)
    lf = L{ii};
    if ~isempty(lf)
        
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        if d(3) < d(2)*1e-10
            jj = 2; % The minumum singular value is degenerate
        else
            jj =3;
        end
        lfo = lf*v(:,jj);
%         w = Cr\lfo / (lfo'/Cr*lfo) ;
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo );
        W{ii} = w;
        VE{ii} = data_filt*w;
        
    end
end

% VE = cell2mat(VE);

ve = data;
ve.trial{1} = cell2mat(VE)';
ve.label = cell(size(ve.trial{1},1),1) ;
for ii = 1:length(ve.label)
    ve.label{ii} = ['ve',num2str(ii)];
end
cfg =[];
cfg.method = 'pca';
cfg.channel = 'all';
comp_pca = ft_componentanalysis(cfg, ve);
score = comp_pca.trial{1}';
compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
icomp = nnz(compvar95) ;
fprintf('%d components for 95perc. of data variance\n',icomp)
clear score

% 
% mri = ft_read_mri('result/PseudoT_win-lose_1s_1-4Hz.nii');
% pseudoT = ft_read_mri('result/PseudoT_win-lose_1s_1-4Hz.nii');
% pseudoT.transform*[ 139 117 79 1]';
% coronal, axial, sagittal
%% Read events

f =data.fsample;

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

% pix = ft_read_data(data_name,'chanindx','UADC016');

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'no';
cfg.channel = 'UADC016';
cfg.demean = 'yes';
pix = ft_preprocessing(cfg);

pixm  = cell2mat(pix.trial)';
pix_white = pixm>2.5;
pix_sample = find(diff(pix_white)==1);



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


win_sample = cell(size(event));
lose_sample = cell(size(event));
ratemood_sample = cell(size(event));
slider_sample = cell(size(event));
gamble_sample = cell(size(event));
certain_sample = cell(size(event));
answer_sample = cell(size(event));
rest_sample = cell(size(event));

for tt = 1:length(event)
    
    if strcmp(event(tt).type, 'UPPT001' ) % Slider appear 3s after        
        switch event(tt).value
            case 1
                answer_sample{tt} = event(tt).sample;
            case 2
                certain_sample{tt} = event(tt).sample;                
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

sample_mismatch = [];

[~,~,m] = find(mood_match.sample);
sample_mismatch = cat(2,sample_mismatch ,m - slider_sample');
[~,~,m] = find(answer_match.sample);
sample_mismatch = cat(2,sample_mismatch ,m - answer_sample');

[~,~,m] = find(choice_match.sample);
m(31)= [];
sample_mismatch = cat(2,sample_mismatch ,m - sort([certain_sample;gamble_sample])');

m = choice_match.sample(logical(choice_match.gamble));
sample_mismatch = cat(2,sample_mismatch ,m - sort(gamble_sample)');

m = outcome_match.sample(logical(abs(outcome_match.win)));
sample_mismatch = cat(2,sample_mismatch ,m - sort([win_sample; lose_sample]'));

%% PLOT
n =4;
comp_topo = zeros(1,length(VE));
comp_topo(grid.inside) = comp_pca.topo(:,n)';
cfg = [];
cfg.parameter = 'pow';
sourceTstat = struct;
sourceTstat.dim = grid.dim;
sourceTstat.inside = grid.inside;
sourceTstat.pos = grid.pos;
sourceTstat.method = 'average';
sourceTstat.avg.pow =  comp_topo;

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


wind_samples =  -data.fsample: 4*data.fsample; % 1s windo
new_trial = find(event_time == 8);
trial_sample = cell(1,length(new_trial)-1);
a = zeros(length(new_trial)-1, length(wind_samples));
for tt =1:length(new_trial)-1
trial_sample{tt} =new_trial(tt) + wind_samples;
a(tt,:) = comp_pca.trial{1}(n,trial_sample{tt});
end
figure; subplot(2,1,1); 
plot(linspace(-1,4,size(a,2)),mean(a,1))


wind_samples =  -2*data.fsample: 2*data.fsample; % 1s windo
new_trial = find(event_time == 70 | event_time == 75 | event_time == 65 );
trial_sample = cell(1,length(new_trial)-1);
a = zeros(length(new_trial)-1, length(wind_samples));
for tt =1:length(new_trial)-1
trial_sample{tt} =new_trial(tt) + wind_samples;
a(tt,:) = comp_pca.trial{1}(n,trial_sample{tt});
end

subplot(2,1,2)
plot(linspace(-2,2,size(a,2)),mean(a,1))

wind_samples =  -2*data.fsample: 2*data.fsample; % 1s windo
new_trial = find(event_time == 55 );
trial_sample = cell(1,length(new_trial)-1);
a = zeros(length(new_trial)-1, length(wind_samples));
for tt =1:length(new_trial)-1
trial_sample{tt} =new_trial(tt) + wind_samples;
a(tt,:) = comp_pca.trial{1}(n,trial_sample{tt});
end

hold on
plot(linspace(-2,2,size(a,2)),mean(a,1))

% %% ICA method
% 
% cfg =[];
% cfg.method = 'fastica';
% cfg.channel = 'all';
% cfg.fastica.numOfIC = icomp;
% comp = ft_componentanalysis(cfg, ve);
% %%
% 
% 
% n = 28;
% comp_topo = zeros(1,length(VE));
% comp_topo(grid.inside) = comp.topo(:,n)';
% cfg = [];
% cfg.parameter = 'pow';
% sourceTstat = struct;
% sourceTstat.dim = grid.dim;
% sourceTstat.inside = grid.inside;
% sourceTstat.pos = grid.pos;
% sourceTstat.method = 'average';
% sourceTstat.avg.pow =  comp_topo;
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
% 
% 
% wind_samples =  -data.fsample: 4*data.fsample; % 1s windo
% new_trial = find(event_time == 8);
% trial_sample = cell(1,length(new_trial)-1);
% a = zeros(length(new_trial)-1, length(wind_samples));
% for tt =1:length(new_trial)-1
% trial_sample{tt} =new_trial(tt) + wind_samples;
% a(tt,:) = comp.trial{1}(n,trial_sample{tt});
% end
% figure; subplot(2,1,1); 
% plot(linspace(-1,4,size(a,2)),mean(a,1))
% 
% 
% wind_samples =  -2*data.fsample: 2*data.fsample; % 1s windo
% new_trial = find(event_time == 70 | event_time == 75 | event_time == 65 );
% trial_sample = cell(1,length(new_trial)-1);
% a = zeros(length(new_trial)-1, length(wind_samples));
% for tt =1:length(new_trial)-1
% trial_sample{tt} =new_trial(tt) + wind_samples;
% a(tt,:) = comp.trial{1}(n,trial_sample{tt});
% end
% 
% subplot(2,1,2)
% plot(linspace(-2,2,size(a,2)),mean(a,1))
% 
% wind_samples =  -2*data.fsample: 2*data.fsample; % 1s windo
% new_trial = find(event_time == 55 );
% trial_sample = cell(1,length(new_trial)-1);
% a = zeros(length(new_trial)-1, length(wind_samples));
% for tt =1:length(new_trial)-1
% trial_sample{tt} =new_trial(tt) + wind_samples;
% a(tt,:) = comp.trial{1}(n,trial_sample{tt});
% end
% 
% hold on
% plot(linspace(-2,2,size(a,2)),mean(a,1))
% %%
% 
% event = ft_read_event(data_name);
% event_time = zeros(data.sampleinfo(2),1);
% % Try localizing button presses in all frequencies to check it is
% % working correctly!!
% for tt = 1:length(event)
%     if strcmp(event(tt).type, 'UPPT001' )
%         event_time(event(tt).sample) = event(tt).value;       
%         
%     end
% end 
% 
% figure; plot(comp.trial{1}(1,:))
% hold on
% plot(event_time*1e-12)