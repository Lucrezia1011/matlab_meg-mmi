function mmi_LTA_subjects_permute(meg_data_name,latent_vars_name,n,npoints,N,fit_parameter,outpathn,datapath)


nn = str2double(n);
npoints = str2double(npoints);

% datapath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/';

meg_data_name = cat(2,datapath,meg_data_name);
latent_vars_name = cat(2,datapath,latent_vars_name);

opts = detectImportOptions(latent_vars_name);
X = readtable(latent_vars_name,opts);

megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 size(X,1)-1])';

w = sqrt(eval(sprintf('X.%s',N)));

Apos = cell(1,N);
Aneg = cell(1,N);
rng('shuffle')

GLM = table;
for z = 1:npoints
    X.MEG = megdata(:,z);
    
    glm = fitglm(X,sprintf('MEG ~ %s',fit_parameter),'Weights',w);  
    
         
    glm_res = glm.Coefficients(2,:);
    glm_res.Row = [];
    glm_res.ROI = nn;
    glm_res.index  = z;
%     eval(['lme_res.',glm_names{ii},' = lme.Coefficients.',glm_names{ii},'(2);'])
  
    GLM = cat(1,GLM,glm_res);
    
end

outpath = [datapath,outpathn,'/glmmodel_',fit_parameter,'/'];
if ~exist(outpath,'dir')
    mkdir(outpath) 
end

outname = [outpath,'ROI_',n,'_permute.txt'];
writetable(GLM,outname)
dlmwrite(outname,[Apos;Aneg]','-append')

end

