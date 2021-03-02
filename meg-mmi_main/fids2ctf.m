function mriC = fids2ctf(mri_name,fids_name,plotOpt)
% Lucrezia Liuzzi 2020/06/29
% Transform nifti mri into ctf coordinates
% mriC = fids2ctf(mri_name,fids_name,plotOpt)
%
% mri_name :  file name of MRI in nift format
% fids_name : fiducial position .tag file
% plotOpt  : =1 for plotting

mri = ft_read_mri(mri_name,'dataformat','nifti');
fileID = fopen(fids_name,'r');
fids_char = fscanf(fileID,'%c');
fclose(fileID);

fids_inds = zeros(3,4);
fids_inds(:,4) = 1;
% 4 lines of 66 characters each
for iF = 1:3
    fids_inds(iF,1) = str2double(fids_char(66*iF+(18:28))); 
    fids_inds(iF,2) = str2double(fids_char(66*iF+(30:40)));
    fids_inds(iF,3) = str2double(fids_char(66*iF+(42:52)));
end

% transformation matrix
T = mri.transform; % = mri.hdr.vox2ras1, shifted by one voxel 
T(1:2,4) = -T(1:2,4); % adjust transformation matrix from fsl to afni (RAI)
T(1,1) = -1; T(2,2) = -1;
% MRI voxels
fid_vox = round(inv(T)*fids_inds')';

% Check with plots ---------------------------------------------------
if plotOpt == 1
    fids_labels = {'Nasion';'Left Ear';'Right Ear'};
    % figure;
    clf
    for iF = 1:3
        subplot(3,3,(iF-1)*3+1)
        imagesc(rot90(squeeze(mri.anatomy(fid_vox(iF,1),:,:)))); axis equal;
        set(gca,'Xdir','reverse')
        title(fids_labels{iF})
        subplot(3,3,(iF-1)*3+2)
        imagesc(rot90(squeeze(mri.anatomy(:,fid_vox(iF,2),:)))); axis equal
        subplot(3,3,(iF-1)*3+3)
        imagesc(fliplr(mri.anatomy(:,:,fid_vox(iF,3)))'); axis equal

    end
    colormap gray
end
% --------------------------------------------------------------------


cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fid_vox(1,1:3); %position of nasion
cfg.fiducial.lpa    = fid_vox(2,1:3); %position of LPA
cfg.fiducial.rpa    = fid_vox(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mriC = ft_volumerealign(cfg,mri);