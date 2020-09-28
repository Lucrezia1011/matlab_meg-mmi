clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

sub = '24071';
freqband  = [13 30];
mu = 0.003;


data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline

if ~exist(data_name,'dir')
    data_name = [sub,'MMI_mmi3_proc1.ds'];
end

if strcmp(sub, '24201') || strcmp(sub, '22695')
    data_name = [sub,'MMI_mmi3_proc2.ds'];
end

cd([data_path,'results'])
mriname = 'ft_coreg_anat.nii';
mri = ft_read_mri(mriname);
mri.coordsys = 'ctf';

%% Segment MRI
cfg = [];
cfg.output  = 'brain';
segmentmri = ft_volumesegment(cfg,mri);

%% Head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentmri);
sens = ft_read_sens(data_name,'senstype','meg');

%% Redefine trials

bv_match = match_triggers_fc(data_name);
slider_match = bv_match.slider;

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = {'UADC005';'UADC006';'UADC007'}; % Time channel!
buttons = ft_preprocessing(cfg);

buttonsd = diff(buttons.trial{1}');
buttonpress = buttonsd>1.5;
[samples,~] = find(buttonpress);


cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);
f = data.fsample;

load([data_path,'ICA_artifacts.mat']);
cfg           = [];
cfg.component = 1:length(comps.label);
data_clean    = ft_rejectcomponent(cfg, comps,data);

% samples = slider_match.sample;
% samples(samples==0) = [];

winds =[-0.5 1];
time =  linspace(winds(1),winds(2),diff(winds)*f+1);
timepoints = winds(1)*f : winds(2)*f;

data_trials = data_clean;

for ii = 1:length(samples)
    data_trials.time{ii} = time;
    data_trials.trial{ii} = data_clean.trial{1}(:,samples(ii)+timepoints);
    data_trials.sampleinfo(ii,:) = samples(ii)+[winds(1)*f, winds(2)*f];
end

%% Find Peak locations

aal_path = '/data/liuzzil2/AAL3/';
load([aal_path,'ROI_MNI_V6_1mm_List.mat'])
aal = ft_read_mri([aal_path,'AAL3_1mm.nii']);

cfg = [];
% cfg.parameter = 'anatomy';
cfg.template = [data_path,'results/ft_coreg_anat.nii'];
% mni =spm_vol([aal_path,'T1.nii']);
mni = ft_read_mri([aal_path,'T1.nii']);
mni.coordsys = 'mni';
% mni.anatomy = mnitemp.anatomy;
[sourcenorm] = ft_volumenormalise(cfg,mni);


cfg = [];
mrinorm = ft_volumenormalise(cfg,mri);

[warped] = ft_warp_apply(mrinorm.cfg.final, input, method, tol)

unix('3dWarp -prefix ft_coreg_anat_deobliq -deoblique ft_coreg_anat.nii')
unix('3dAllineate -base ft_coreg_anat_deobliq+orig -input T1.nii -output T1shift')

aal_path = '/data/liuzzil2/AAL3';
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
