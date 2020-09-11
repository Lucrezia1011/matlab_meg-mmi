
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit

out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';
meg_data_name = sprintf('%s/powergrid_25-40Hz.txt',out_path);
% meg_data_name = sprintf('%s/powersens_25-40Hz.txt',out_path);

latent_vars_name = [out_path,'/latent_vars.csv'];

%%
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';


datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_25-40Hz');

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/mni_grid.txt');
Te = zeros(size(gridall));

param_list = cell(1,15);
for nn = 1:14
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end


cd([datapath,'/lme_E'])

X = zeros(nnz(gridall),1);
pV =  zeros(nnz(gridall),1);
for nn = 1:14
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
end
Te(gridall==1) = X;
%Step 2 is met: the causal variable is correlated with the mediator
% Bonferroni
zz = find(pV<0.05/nnz(pV));


nullnames = dir('grid_permute');
nullnames(1:2)=[];
clusternull = zeros(length(nullnames),nnz(gridall));
for n = 1:length(nullnames)
    clusternull2 = dlmread(['grid_permute/',nullnames(n).name]);
    clusternull(n,:) = clusternull2;
end

E = 0.5; %0.1  % try and change the parameters
H = 2; %2
dh = 0.1; 
M = zeros(1,length(nullnames));
for n = 1:length(nullnames)
    clusternull2 = zeros(size(gridall));
    clusternull2(gridall==1) = clusternull(n,:);
    img = reshape(clusternull2,sourcemodel.dim);
    
    tfce= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
    M(n) = max(tfce(:));
    
end
M = sort(M,'descend');
thresh = M(0.05*size(clusternull,1)); % 0.01
    
img = reshape(Te,sourcemodel.dim);

tfce= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
tfce = tfce(gridall==1);

zz = find(tfce>thresh);

%%

opts = detectImportOptions(latent_vars_name);
ltv = readtable(latent_vars_name,opts);

% fid=fopen(meg_data_name);
% g = textscan(fid,'%s','delimiter','\n');
% fclose(fid);
% ntot = length(g{1});
% clear g
megdata = dlmread(meg_data_name)';

% Includes intercep random effect for subject
% Uncorrelated random effect for intercept and mood grouped by subject,
% plus random effect for intercept given trial

% lme_xm = sprintf('MEG ~ E +(E|subject) + (1|trial) + (1|recording)');
% lme_xy = sprintf('mood ~ E + (E|subject) + (1|trial)');
% lme_xmy = sprintf('mood ~ MEG + E +(E|subject) + (1|trial) + (-1 + MEG|recording)');

lme_xm = sprintf('E ~ MEG + (1|subject) + (-1 + MEG|trial) + (-1 + MEG|recording)');
lme_xy = sprintf('mood ~ MEG + (1|subject) + (-1 + MEG|trial) + (-1 + MEG|recording)');
lme_xmy = sprintf('mood ~ MEG + E + (E|subject) + (-1 + MEG|trial) + (-1 + MEG|recording)');

% X.mood = zscore(X.mood); % cannot z-score all together!
% X.E = zscore(X.E);
subs = unique(ltv.subject)';
for s = subs
    ltv.mood(ltv.subject==s) = zscore(ltv.mood(ltv.subject==s));
    ltv.E(ltv.subject==s) = zscore(ltv.E(ltv.subject==s));
end
%%
M = cell(1,length(zz));
parfor z = 1:length(zz)
    X = ltv;
    X.MEG = (megdata(:,zz(z)));
    for r = 1:max(X.recording)
        X.MEG(X.recording==r) = zscore(X.MEG(X.recording==r));
    end
    
    lme1 = fitlme(X,lme_xy); %Expectation effect on mood
    lme2 = fitlme(X,lme_xm); %Expectation effect on MEG (i.e. the mediator)
    lme3 = fitlme(X,lme_xmy); %Expectation and mediator effects on mood
    
    LME = [];
    % Step 1
    LME.xy_tStat = lme1.Coefficients.tStat(2); % Equivalent of Estimate/SE
    LME.xy_pValue = lme1.Coefficients.pValue(2);
    
    % Step 2
    LME.xm_tStat = lme2.Coefficients.tStat(2);
    LME.xm_pValue =lme2.Coefficients.pValue(2);
    
    % Step 3
%     LME.xym_Name = lme3.Coefficients.Name;
    LME.xym_tStat = lme3.Coefficients.tStat(3);
    LME.xym_pValue = lme3.Coefficients.pValue(3);
    
    LME.c = lme1.Coefficients.Estimate(2);
    LME.a = lme2.Coefficients.Estimate(2);
    LME.c1 = lme3.Coefficients.Estimate( strcmp(lme3.Coefficients.Name,'MEG'));
    LME.b = lme3.Coefficients.Estimate( strcmp(lme3.Coefficients.Name,'E'));
   
    
    % The value of the mediated or indirect effect estimated by taking the 
    % difference in the coefficients, ĉ − ĉ′, corresponds to the reduction 
    % in the independent variable effect on the dependent variable when 
    % adjusted for the mediator. To test for significance, the difference 
    % is then divided by the standard error of the difference and the ratio 
    % is compared to a standard normal distribution.
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819368/
    
    % Step 4
    LME.mediation = (LME.c-LME.c1)/sqrt(lme1.Coefficients.SE(2)^2 + lme3.Coefficients.SE( strcmp(lme3.Coefficients.Name,'MEG'))^2);
    LME.ab = LME.a*LME.b;
   
    % Sobel Test:  SE_{ab} = b^2*SEb^2 + a^2*SEa^2
    LME.SE_ab = sqrt( LME.a^2 * lme2.Coefficients.SE(2)^2 + ...
        LME.b^2 * lme3.Coefficients.SE( strcmp(lme3.Coefficients.Name,'E'))^2 );
    
   
%   The amount of mediation is called the indirect effect.   Note that the
%   total effect = direct effect + indirect effect; 
%   or using symbols:  c = c' + ab
%  the two are only approximately equal for multilevel models, logistic 
%  analysis and structural equation modeling with latent variables.  
%  For such models, it is probably inadvisable to compute c from Step 1, 
%  but rather c or the total effect should be inferred to be c' + ab and 
%  not directly computed.


% Joint Significance of Paths a and b: Test that a and b are not zero.
% i.e. xm_tStat & lme3.Coefficients.tStat( strcmp(lme3.Coefficients.Name,'E'))
    LME.jointSig = LME.xm_pValue<0.05 & ...
        lme3.Coefficients.pValue( strcmp(lme3.Coefficients.Name,'E'))<0.05;
        
    M{z} = LME;
    clc
    fprintf('Done %.1f\n',z/length(zz)*100)
end

save('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_25-40Hz/lme_E/mediation_analysis.mat','M','zz');

%Most contemporary analysts believe that the essential steps in establishing 
% mediation are Steps 2 and 3. Certainly, Step 4 does not have to be met 
% unless the expectation is for complete mediation.  In the opinion of most 
% though not all analysts, Step 1 is not required.

% Monte Carlo http://quantpsy.org/medmc/medmc.htm

% indirect effect
ab = zeros(length(zz),2);
for z = 1:length(zz)
    % Z-test > 1.96 in absolute value is significant at the .05 leve
    ab(z,1) = M{z}.ab/M{z}.SE_ab;
    ab(z,2) = M{z}.jointSig;
end

AB = zeros(size(pV));
AB(zz) = ab(:,1);
img = zeros(size(Te));
img(gridall==1) = AB;

mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');

sourceant.pow  = img;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';

crang = [1.96 4];
% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'slice'; %'ortho'
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
cfg.method        = 'ortho'; %'ortho'
ft_sourceplot(cfg, sourceout_Int);


