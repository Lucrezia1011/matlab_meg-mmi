clear all
close all
clc
meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

data_list = [];

zz= 0;
% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
Nlist([10,24,26,49,53]) = [];
for sn = Nlist
        
    sdan = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;
            data_list{zz} = data_name;
        end
    end
    
   
end

%%
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
% 
% % roiopt = 'c';
roiopt = 'grid';
% % tfsopt = 'pwelch';
% tfsopt = 'm';
% 
gridres= 5;


freqband = [13 30]; 
windLength = 0.5;  
windWait = 0.6;
mu =4;
plotopt = 'anat';
for ii = 58:59% 14:19%length(data_list)
    mmiButtonpressBF(data_list{ii},freqband,mu,gridres,windLength,windWait,plotopt)
end


%% Read results
Tstat_all = cell(1,length(data_list));

for ii =1:length(data_list)
    data_name = data_list{ii};
    sub = data_name(5:9);
    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
    save_name = sprintf('%sbuttonpressBF_%.0f-%.0fHz_a%.1fs_c%.1fs.mat',...
        processing_folder,freqband(1),freqband(2),windLength,windWait);
    load(save_name);
    Tstat_all{ii} = Tstat;
end
%% Plot with MNI
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');

for ii = 13%1:12
    
Tstat = Tstat_all{ii};
sourceant = [];
sourceant.pow = Tstat;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';


crang = [];
cfg = [];
cfg.method        = 'slice'; %'slice'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
ft_sourceplot(cfg, sourceout_Int);
set(gcf,'name',data_list{ii})
end


%% Plot in anatomical space
ii = 13;
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
% mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');

data_name = data_list{ii};
sub = data_name(5:9);
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)
mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

Tstat = Tstat_all{ii};


template_grid = sourcemodel;
clear sourcemodel
% sourcemodel based on 5mm grid MNI brain
cfg = [];
cfg.mri = mri;
cfg.warpmni = 'yes';
cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
cfg.unit      = 'm';
cfg.nonlinear = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);


sourceant = [];
sourceant.pow = Tstat;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'ctf';

crang = [min(sourceout_Int.pow(:)) min(sourceout_Int.pow(:))/2];
% crang = [-.3 -.1];
cfg = [];
cfg.method        = 'ortho'; %'slice'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
ft_sourceplot(cfg, sourceout_Int);
