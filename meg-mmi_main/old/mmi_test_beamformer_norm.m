clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

freqband = [40 90];
ii = 0;
mri_norm = cell(1,11);
source_norm = cell(1,11);
zname_nii_cat = [];
task_name = 'Buttonpress_0.1s';
% task_name = 'Mooddiff_quartile';


for s = [1,3,4,6,7,8,9,11,14,15,16]% [6,8,9,14] %[3,6,7,8,9,11,14] %[1,4,15,16]  %[1,3,4,6,7,8,9,11,14,15,16] %Subjects showing enough variation in mood
    sub = subn(s,:);
    
   
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    if ~exist(data_name,'dir')
        data_name = [sub,'MMI_mmi3_proc1.ds'];
    end
    
    if strcmp(sub, '24201') || strcmp(sub, '22695')
        data_name = [sub,'MMI_mmi3_proc2.ds'];
    end
    
    
    mri_name = [sub,'_anat+orig.BRIK'];

    %% Template transform
    cd results
    zname =  sprintf('ft_coreg_anat_norm');
    zname_nii = [zname,'.nii'];
    mrinorm = ft_read_mri(zname_nii,'dataformat','nifti');
    
    mu = 0.001;
%     zname =  sprintf('Mooddiff_3s_PseudoT_%d-%dHz_norm_mu%.0s',freqband(1),freqband(2),mu);
    zname =  sprintf([task_name,'_PseudoT_%d-%dHz_norm_mu%.0s'],freqband(1),freqband(2),mu);
    zname_nii = [zname,'.nii'];
    sourcenorm = ft_read_mri(zname_nii,'dataformat','nifti');
    
    ii = ii+1;
    mri_norm{ii} = mrinorm;
    source_norm{ii} = sourcenorm;

    zname_nii_cat = cat(2,zname_nii_cat,' ',[data_path,'results/',zname_nii]);
end


%%
anat  =zeros([sourcenorm.dim,ii]);
Tstat = zeros([sourcenorm.dim,ii]);
for jj = 1:ii
    anat(:,:,:,jj) = mri_norm{jj}.anatomy;
    Tstat(:,:,:,jj) = source_norm{jj}.anatomy;
end
Tstatm = mean(Tstat,4);
anatm = mean(anat,4);
% 
% 
% figure; imagesc(Tstatm(:,:,50)); caxis([-.15 .15])

%% Write average

cd /data/MBDU/MEG_MMI3/results

zname_pre = sprintf([task_name,'_PseudoT_%d-%dHz_mu%.0s_%dsingles.nii'],freqband(1),freqband(2),mu,ii);
unix(['3dTcat -tpattern seqplus -prefix ',zname_pre,' ',zname_nii_cat]);



mri_norm_av = mrinorm;
mri_norm_av.anatomy = anatm;

Tstat_norm_av = sourcenorm;
Tstat_norm_av.anatomy = Tstatm;


% ft_write_mri('ft_coreg_anat_norm_11av.nii',mri_norm_av,'dataformat','nifti')
% ft_write_mri(sprintf('Mooddiff_3s_PseudoT_%d-%dHz_mu%.0s_%dav.nii',freqband(1),freqband(2),mu,ii),...
%     Tstat_norm_av,'dataformat','nifti')
ft_write_mri(sprintf([task_name,'_PseudoT_%d-%dHz_mu%.0s_%dav.nii'],freqband(1),freqband(2),mu,ii),...
    Tstat_norm_av,'dataformat','nifti')



% 4 HV average
s = [4,6,7,9]; %[6,8,9,14];
Tstat_norm_av = sourcenorm;
Tstat_norm_av.anatomy = mean(Tstat(:,:,:,s),4);

ft_write_mri(sprintf([task_name,'_PseudoT_%d-%dHz_mu%.0s_%davHV.nii'],freqband(1),freqband(2),mu,4),...
    Tstat_norm_av,'dataformat','nifti')

% 4 MDD average
s = [1,3,10,11];%[1,4,15,16];
Tstat_norm_av = sourcenorm;
Tstat_norm_av.anatomy = mean(Tstat(:,:,:,s),4);

ft_write_mri(sprintf([task_name,'_PseudoT_%d-%dHz_mu%.0s_%davMDD.nii'],freqband(1),freqband(2),mu,4),...
    Tstat_norm_av,'dataformat','nifti')

