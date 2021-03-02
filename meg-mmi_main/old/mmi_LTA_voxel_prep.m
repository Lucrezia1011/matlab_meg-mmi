function mmi_LTA_voxel_prep(data_name,twind,evokedopt,inducedopt)
% [1:7,11,14:16]
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
addpath ~/matlab

%% Co-register MRI from fiducial positions

% LTA model latent variables:
% EC: Expectation of certain value
% EG: Expectation during gabling
% Ediff: Drift rate
% LTA: Long term average with gamma:   1/t * sum_i=1 ^t(V(i)^gamma),   cumsum(LTA.OutcomeAmount^gamma)./(1:ntrials)'
% V_i^gamma = outcome of trial i
% new_p = subjective winning probability
% RPE = Reward prediction error
% LTA_sum  = sum(LTA)
% RPE_sum = sum(RPE)
% log_like
% mood_log_like

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)


bv_names = dir('/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/');
for ii = 1:length(bv_names)
    if strcmp(bv_names(ii).name,['LTA_Gamma_Latent_Variables_3_Blocks_MEG_',sub,'.csv'])
        bv_name = ['/data/MBDU/MEG_MMI3/data/LTA_Gamma_181219/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
ltv = readtable(bv_name,opts); % bahavioral data



%% Co-register MRI

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if ~exist([processing_folder,'/leadfields_5mm_MNI.mat'],'file')

    %% AAL atlas
    gridres = 5; % resolution of beamformer grid in mm

    % Load fieldtrip 10mm MNI grid
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    template_grid = sourcemodel;

    % atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
    % atlas = ft_convert_units(atlas,sourcemodel.unit);
    % 
    % cfg = [];
    % cfg.interpmethod = 'nearest';
    % cfg.parameter = 'tissue';
    % sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

    clear sourcemodel

%% Sourcemodel warp MNI grid

    % sourcemodel based on 5mm grid MNI brain
    cfg = [];
    cfg.mri = mri;
    cfg.warpmni = 'yes';
    cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
    cfg.unit      = 'm';
    cfg.nonlinear = 'yes';
    sourcemodel = ft_prepare_sourcemodel(cfg);

%% Find location of AAL ROIs
% R = length(sourcemodelAAL.tissuelabel);
% locs = zeros(R,3);
% locsAAL = zeros(R,3);
% for ii = 1:R
%     ind = find(sourcemodelAAL.tissue == ii);
%     locs(ii,:) = mean(sourcemodel.pos(ind,:));
%     
%     locsAAL(ii,:) = mean(template_grid.pos(ind,:));
% end


%% Calculate lead fields
    cfg                 = [];
    cfg.grad            = sens;
    cfg.headmodel       = vol;
    cfg.reducerank      = 2;
    cfg.channel         = {'MEG'};
    cfg.sourcemodel.pos = sourcemodel.pos;
    cfg.sourcemodel.unit   = 'm';
    cfg.siunits         = true;
    cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
    [grid] = ft_prepare_leadfield(cfg);
    save([processing_folder,'/leadfields_5mm_MNI.mat'],'grid');
else
    load([processing_folder,'/leadfields_5mm_MNI.mat']);
end
%% Clean data with ICA

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
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data          = ft_rejectcomponent(cfg, comps,data);

%%
filt_order = []; % default
%     data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[0.5 150],filt_order,'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample,35,filt_order,'but');

data.trial{1} = data_filt;
clear data_filt


%% Read events

bv_match = match_triggers_fc(data_name);

outcome_match  = bv_match.outcome;
cue_match = bv_match.answer;
choice_match = bv_match.choice;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;


%% Mood

indbm = blockmood_match.sample~=0;
indm = mood_match.sample~=0;

[x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
Fmood = griddedInterpolant(x,v(ind),'pchip');

xi = outcome_match.sample(outcome_match.win~=0);
moodi = Fmood(xi); % interpolated mood timecourse
 

%% Beamfomer
icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

% Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value

L = grid.leadfield;

% VE(1:size(L,2)) = {0};
% W(1:size(L,2)) = {0};


datave = struct;
datave.trial{1} = data.trial{1}(1,:);
datave.label = {'ve'};
datave.time = data.time;
datave.fsample = f;
datave.sampleinfo = data.sampleinfo;

vox = find(grid.inside);
VE_cue = cell(1,nnz(vox));
VE_choice = cell(1,nnz(vox));
VE_out = cell(1,nnz(vox));

for ii = 1:length(vox)
    
    lf = L{vox(ii)}; % Unit 1Am
    
    
    % %  G O'Neill method, equivalent to ft
    [v,d] = svd(lf'/Cr*lf);
    d = diag(d);
    jj = 2;
    
    lfo = lf*v(:,jj); % Lead field with selected orientation
    
    w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
    
    %         VE =  w'*data.trial{1};
    %         W{ii} = w;
    %         VE{ii}  = w'*data.trial{1};
    datave1 = datave;
    datave1.trial{1} =  w'*data.trial{1};
    
    [dataout,~]= define_trials(outcome_match.sample(outcome_match.win~=0), datave1, tasktime, [-.2,1],0);
    % Anticipatory window lasts 4s, unlikely to see evoked responses:
    % try oscillatory!
    %         cfg = [];
    %         cfg.resamplefs = 300; % Downsample to 200Hz for ease of memory
    %         cfg.feedback = 'no';
    %         dataout = ft_resampledata(cfg, dataout);
    
    dataout = cell2mat(dataout.trial');
    VE_out{ii} = resample(dataout',300,f);
    
    [datacue,~]= define_trials(cue_match.sample(cue_match.sample~=0), datave1, tasktime, [-.2,1],0);
    % Anticipatory window lasts 4s
    datacue = cell2mat(datacue.trial');
    VE_cue{ii} = resample(datacue',300,f);
    
    [datachoice,~]= define_trials(choice_match.sample(choice_match.sample~=0), datave1, tasktime, [-.2,1],0);
    datachoice = cell2mat(datachoice.trial');
    VE_choice{ii} = resample(datachoice',300,f);
    %
    %         figure; subplot(2,3,1)
    %         plot(datacue.time{1},mean(cell2mat(datacue.trial'),1))
    %         title('cue')
    %         subplot(2,3,2)
    %         plot(datachoice.time{1},mean(cell2mat(datachoice.trial'),1))
    %         title('choice')
    %         subplot(2,3,3)
    %         plot(dataout.time{1},mean(cell2mat(dataout.trial'),1))
    %         title('feedback')
    
    % TFS
    %
    %         [dataout,ttdel]= define_trials(outcome_match.sample(outcome_match.win~=0), datave, tasktime, [-4,4]);
    %         % Anticipatory window lasts 4s, unlikely to see evoked responses:
    %         % try oscillatory!
    %
    %         [datacue,ttdel]= define_trials(cue_match.sample(cue_match.sample~=0), datave, tasktime, [-4,4]);
    %
    %         [datachoice,ttdel]= define_trials(choice_match.sample(choice_match.sample~=0), datave, tasktime, [-4,4]);
    %
    %
    %
    %         cfg  =[];
    %         cfg.taper = 'hanning';
    %         cfg.output = 'pow';
    %         cfg.pad = 'nextpow2';
    %         cfg.method = 'mtmconvol';
    %         cfg.foi    = 0.5:0.5:35;
    %         cfg.t_ftimwin    = 7./cfg.foi;
    %         cfg.keeptrials = 'yes';
    %         cfg.toi    = -4:0.05:4;
    %
    %         tfs_out   = ft_freqanalysis(cfg,dataout);
    %         tfs_cue   = ft_freqanalysis(cfg,datacue);
    %         tfs_choice   = ft_freqanalysis(cfg,datachoice);
    %
    %         cfgf = [];
    %         cfgf.baseline     = [-4 4];
    %         cfgf.baselinetype = 'absolute';
    %         cfgf.maskstyle   = 'saturation';
    %         cfgf.interactive  = 'no';
    %         cfgf.zlim         = [-1 1]*3e-27;
    %         cfgf.xlim         = [-2 2];
    % %         cfgf.ylim        = [1 35];
    %         cfgf.trials  = 'all';
    %
    %         subplot(2,3,4)
    %         ft_singleplotTFR(cfgf, tfs_cue);
    %         subplot(2,3,5)
    %         ft_singleplotTFR(cfgf, tfs_choice);
    %         subplot(2,3,6)
    %         ft_singleplotTFR(cfgf, tfs_out);
      
    if mod(ii,1000) == 0
        clc
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
    
end
clc

[~,ttdel_out]= define_trials(outcome_match.sample(outcome_match.win~=0), datave, tasktime, [-.2,1],0);
[~,ttdel_cue]= define_trials(cue_match.sample(cue_match.sample~=0), datave, tasktime, [-.2,1],0);
[~,ttdel_choice]= define_trials(choice_match.sample(choice_match.sample~=0), datave, tasktime, [-.2,1],0);


ltvind_out = outcome_match.bv_index(outcome_match.win~=0) - 12; % indeces start at 13
ltvind_out(ttdel_out) = [];

ltvind_cue = cue_match.bv_index(cue_match.sample~=0) - 12; % indeces start at 13
ltvind_cue(ttdel_cue) = [];

ltvind_choice = choice_match.bv_index(choice_match.sample~=0) - 12; % indeces start at 13
ltvind_choice(ttdel_choice) = [];

save([processing_folder,'/mmi_voxel_prep.mat'],'-v7.3','VE_cue','VE_choice','VE_out');
%     fprintf('Beamformer finished\n' )

% VE = cell2mat(VE');

%% Align evoked responses
VE_align
findNearestNeighbors
for ii = 1:length(L)
    d = sqrt(sum((grid.pos-grid.pos(ii,:)).^2,2));
    
    grid.pos = 1 ;% Unit 1Am
    
end

%%

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end


% Try this with all subjects!

S = str2double(sub)*ones(ntrials,1);
% sen = repmat(dataout.label,[ntrials,1]);

ltvcut = [ltv.EC, ltv.EG, ltv.Ediff, ltv.LTA, ltv.new_p, ltv.RPE, ltv.LTA_sum, ltv.RPE_sum];
ltvcut = ltvcut(ltvind,:);
trials = ltv.Var1;
trials = trials(ltvind);

% ltvcut(:,4:11) = ltvcut;
% ltvcut(:,1) = S;
% ltvcut(:,2) = trials;
% ltvcut(:,3) = moodi;
moodi(ttdel) = [];
ltvall = table(S,trials,moodi',ltvcut(:,1),ltvcut(:,2),ltvcut(:,3),...
    ltvcut(:,4),ltvcut(:,5),ltvcut(:,6),ltvcut(:,7),ltvcut(:,8),'VariableNames',...
        {'subject','trial','mood','EC','EG','Ediff','LTA','new_p','RPE','LTA_sum','RPE_sum'});


% Y = cell2mat(dataout.trial);
% Y = reshape(Y,[nchans, size(dataout.trial{1},2),length(dataout.trial)]);

save(save_name ,'Y','ltvall','-append')


% Do this at the end??
%     % for evoked responses!! To deal with sign uncertainty
%     if ~isempty(Yall)
%         vec = corr(mean(Y,3)');
%         Y = Y.*sign(vec(:,69)); % aligns to the left precentral lobule
%         vec = corr(mean(Y,3)',mean(Yall,3)');
%         Y = Y.*sign(vec(69,69));
%     else
%         vec = corr(mean(Y,3)');
%         Y = Y.*sign(vec(:,69));
%     end
%%%%%%%%%%%%%%
end



