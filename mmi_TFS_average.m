clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

% ii = 0;
% mri_norm = cell(1,11);
% source_norm = cell(1,11);
% zname_nii_cat = [];

freqnames = {'delta';'theta';'alpha';'beta';'gamma_low';'gamma_high'};
mrinorm = ft_read_mri('/home/liuzzil2/fieldtrip-20190812/external/spm8/templates/T1.nii');
for k = 1:length(freqnames)
    ss = freqnames{k};
    mrinorm.(ss) = [];
end
for sn = 1:16 % [6,8,9,14] %[3,6,7,8,9,11,14] %[1,4,15,16]  %[1,3,4,6,7,8,9,11,14,15,16] %Subjects showing enough variation in mood
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    if exist('Mood_corr_alpha.nii','file')
    mri = ft_read_mri([sub,'_coreg.nii']);
   
    for k = 1:length(freqnames)
        Timage = ft_read_mri(['Mood_corr_',freqnames{k},'.nii']);
        ss = freqnames{k};
        mri.(ss) = fisher(Timage.anatomy); % Transfor r values to z
    end
    mri.coordsys = 'ctf';
    
    cfg = [];
%     cfg.template = '/home/liuzzil2/MNI152_T1_2009c.nii';
    [sourcenorm] = ft_volumenormalise(cfg, mri);
    
    for k = 1:length(freqnames)
        ss = freqnames{k};
        mrinorm.(ss) = cat(4,mrinorm.(ss),sourcenorm.(ss));
    end
    
    end
%     mrinorm = ft_read_mri('/home/liuzzil2/MNI152_T1_2009c.nii');
%     mrinorm.delta = sourcenorm.delta;
end

sn = 4;
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
mri = ft_read_mri([sub,'_coreg.nii']);
mri.coordsys = 'ctf';
cfg = [];
[sourcenorm] = ft_volumenormalise(cfg, mri);
    
mriplot = sourcenorm;
for k = 1:length(freqnames)
    ss = freqnames{k};
    sn_split = [7,10,12];
    sn_cont = [1:6,8,9,11];
    sn_lowvar = [2,5];
    sn_bigvar = [3,6,8,9,11];
    
    sn_cont_highvar = [1,3,4,6,8,9,11];
    
    mriplot.(ss) = mean(mrinorm.(ss)(:,:,:,sn_cont_highvar),4);  
end


freqnamest = {'delta';'theta';'alpha';'beta';'gamma_{low}';'gamma_{high}'};

for k = [2,4,5] %1:length(freqnames)
    ss = freqnames{k};
%     crang = [-0.4 -0.1];
%     cfg = [];
%     cfg.method        = 'slice';
%     cfg.funparameter = ss;
%     cfg.maskparameter = ss;
%     cfg.funcolormap  = 'auto';
%     cfg.funcolorlim   = crang;
%     cfg.opacitylim = crang;
%     ft_sourceplot(cfg, mriplot);
%     title(freqnamest{k})
    
    crang = [0.1 0.4];
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter = ss;
    cfg.maskparameter = ss;
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    ft_sourceplot(cfg, mriplot);
    title(freqnamest{k})
    
end


%% Check high gamma correlation
sn_cont_highvar = [1,3,4,6,8,9,11];
k = 6;

for sn = sn_cont_highvar
    mriplot.(ss) = mrinorm.(ss)(:,:,:,sn);
 
    ss = freqnames{k};
    crang = [-0.4 -0.1];
    cfg = [];
    cfg.method        = 'slice';
    cfg.funparameter = ss;
    cfg.maskparameter = ss;
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    ft_sourceplot(cfg, mriplot);
    title(['Subject ',num2str(sn),' ',ss])
    
    
end

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

