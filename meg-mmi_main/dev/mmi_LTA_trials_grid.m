function mmi_LTA_trials_grid(meg_data_name,latent_vars_name,n,fit_parameter,outpath)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit

opts = detectImportOptions(latent_vars_name);
X = readtable(latent_vars_name,opts);

% do 336 points at a time
% megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 size(X,1)-1])';
megdata = dlmread(meg_data_name)';
% megdata = reshape(megdata,[60,336,size(megdata,3)] );

% Includes intercep random effect for subject
% Uncorrelated random effect for intercept and mood grouped by subject,
% plus random effect for intercept given trial

% Previously potentially correlated random effect for intercept and mood grouped by subject.
%     lme = fitlme(X,sprintf('MEG ~ %s + (%s - 1 | subject) + (1|trial)',fit_parameter,fit_parameter));
% include recordig random variable for tfs
lme_formula = sprintf('MEG ~ %s +(%s|subject) + (1|trial) + (1|recording)',fit_parameter,fit_parameter);
% lme_formula = sprintf('MEG ~ %s + (%s | subject) + (1|trial)',fit_parameter,fit_parameter);

npoints = size(megdata,2);
% LME = table;
T = zeros(npoints,1);
for z = 1:npoints % 60 time points * 336 voxels
    X.MEG = megdata(:,z);
    
    lme = fitlme(X,lme_formula); %Fixed effects for RPE
    
%     lme_names = lme.Coefficients.Properties.VarNames;
     
%     lme_res = table;
%     lme_res.vox = nn;
%     lme_res.time = mod(z,60);
%     for ii = 2:length(lme_names)
%         lme_res.(lme_names{ii}) = lme.Coefficients.(lme_names{ii})(2);
% %         eval(['lme_res.',lme_names{ii},' = lme.Coefficients.',lme_names{ii},'(2);'])
%     end
    
%     LME = cat(1,LME,lme_res);
    T(z) = lme.Coefficients.tStat(2); % fit_paramter predictor
    
end

% outpath = [data_path,outpathn,'/lme_',fit_parameter,'/'];
% if ~exist(outpath,'dir')
%     mkdir(outpath) 
% end

outname = [outpath,'voxs_',n,'.txt'];
% writetable(LME,outname)
dlmwrite(outname,T);
end