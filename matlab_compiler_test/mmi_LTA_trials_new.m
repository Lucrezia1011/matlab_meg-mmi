function mmi_LTA_trials_new(meg_data_name,latent_vars_name,n,npoints,fit_parameter,outpath)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit
nn = str2double(n);
npoints = str2double(npoints);

% datapath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/';

% meg_data_name = cat(2,data_path,meg_data_name);
% latent_vars_name = cat(2,data_path,latent_vars_name);

opts = detectImportOptions(latent_vars_name);
X = readtable(latent_vars_name,opts);

megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 size(X,1)-1])';

% Includes intercep random effect for subject
% Uncorrelated random effect for intercept and mood grouped by subject,
% plus random effect for intercept given trial

% Previously potentially correlated random effect for intercept and mood grouped by subject.
%     lme = fitlme(X,sprintf('MEG ~ %s + (%s - 1 | subject) + (1|trial)',fit_parameter,fit_parameter));
if npoints == 925
    % include recordig random variable for tfs
    lme_formula = sprintf('MEG ~ %s + (-1 + %s | subject) + (1|recording) + (1|trial)',fit_parameter,fit_parameter); 
else
    lme_formula = sprintf('MEG ~ %s + (%s | subject) + (1|trial)',fit_parameter,fit_parameter); 
end

LME = table;
for z = 1:npoints
    X.MEG = megdata(:,z);
    
    lme = fitlme(X,lme_formula); %Fixed effects for RPE
    
    lme_names = lme.Coefficients.Properties.VarNames;
     
    lme_res = table;
    lme_res.ROI = nn;
    lme_res.index = z;
    for ii = 2:length(lme_names)
        lme_res.(lme_names{ii}) = lme.Coefficients.(lme_names{ii})(2);
%         eval(['lme_res.',lme_names{ii},' = lme.Coefficients.',lme_names{ii},'(2);'])
    end
    LME = cat(1,LME,lme_res);
    
end

% outpath = [data_path,outpathn,'/lme_',fit_parameter,'/'];
% if ~exist(outpath,'dir')
%     mkdir(outpath) 
% end

outname = [outpath,'ROI_',n,'.csv'];
writetable(LME,outname)

end