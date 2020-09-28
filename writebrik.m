function writebrik(name,mri)
ft_write_mri([name,'.nii'],mri,'dataformat','nifti');
unix(['3dcopy ',name,'.nii ',name])
unix(['3drefit -orient RAI ',name,'+orig'])
unix(['3dWarp -deoblique -prefix ',name,'_ortho ',name,'+orig'])
unix(['3drefit -orient ALI ',name,'_ortho+orig'])

