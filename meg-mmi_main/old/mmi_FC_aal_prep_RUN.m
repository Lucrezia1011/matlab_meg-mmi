clear all
close all
clc

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn = 1:16 %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') ...
                && exist([data_path,name_list(ii).name,'/beamforming/ICA_artifacts.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)   
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

%%
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
addpath ~/ppyll1/matlab/svdandpca

twind = 60;
inducedopt = {'gamma'};%{'theta';'alpha';'beta'}; % cell with freq limitis or empty
roiopt = 'g';
Popt = 0;
FCopt =1;
for ii = 4 %1:length(param_list) 
        mmi_FC_aal_prep(param_list{ii},twind,inducedopt,roiopt,Popt,FCopt)  
%         mmi_FC_aal_fix(param_list{ii},twind,inducedopt,roiopt)  
end

return
%%

ltvMood = [];

YTheta = [];
YBeta = [];
YAlpha = [];
YDelta = [];
YGamma = [];
z = 0;
for sn = 1:16
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') ...
                && exist([data_path,name_list(ii).name,'/beamforming/ICA_artifacts.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)  
        z = z+1;
        data_name = data_names{runs};
        n = str2double(data_name(end-3));
        if isnan(n)
            n = 1;
        end
        
        load(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/FC/',sub,'_',num2str(n)])

        YDelta = cat(2,YDelta,Pow_delta);
        YTheta = cat(2,YTheta,Pow_theta);
        YAlpha = cat(2,YAlpha,Pow_alpha);
        YBeta = cat(2,YBeta,Pow_beta);
        YGamma = cat(2,YGamma,Pow_gamma);
        
        ltvmood.recording = repmat(z,size(ltvmood,1),1);
        if isempty(ltvMood)
            ltvMood = ltvmood;  
        else
            ltvMood(end+(1:size(ltvmood,1)),:) = ltvmood;
        end
                
    end
end

clear ltvmood

%% Linear mixed effects model

% Run linear mixed model with frequency power envelope
nrois = size(YAlpha,1);
% beta_coeff = zeros(2,nrois);
% delta_coeff = zeros(2,nrois);
% alpha_coeff = zeros(2,nrois);
% theta_coeff = zeros(2,nrois);
% gamma_coeff = zeros(2,nrois);

beta_coeff = cell(1,nrois);
delta_coeff = cell(1,nrois);
alpha_coeff = cell(1,nrois);
theta_coeff = cell(1,nrois);
gamma_coeff = cell(1,nrois);

% lme_formula = 'MEG ~ mood +(-1 + mood|subject) + (1|trial) + (1|recording)';
lme_formula = 'MEG ~ trial + (trial|subject) + (1|recording)';

% what if I tried to find a fixed effect of trial?
parfor n = 1:nrois
    X = ltvMood;
    
    X.MEG = YDelta(n,:)';
    lme_delta = fitlme(X,lme_formula); 
    delta_coeff{n} = [lme_delta.Coefficients.tStat(2); ...
        lme_delta.Coefficients.pValue(2)];
    
    X.MEG = YTheta(n,:)';
    lme_theta = fitlme(X,lme_formula); 
    theta_coeff{n} = [lme_theta.Coefficients.tStat(2); ...
        lme_theta.Coefficients.pValue(2)];
    
    X.MEG = YAlpha(n,:)';
    lme_alpha = fitlme(X,lme_formula); 
    alpha_coeff{n} = [lme_alpha.Coefficients.tStat(2); ...
        lme_alpha.Coefficients.pValue(2)];
    
    X.MEG = YBeta(n,:)';
    lme_beta = fitlme(X,lme_formula); 
    beta_coeff{n} = [lme_beta.Coefficients.tStat(2); ...
        lme_beta.Coefficients.pValue(2)];
    
    X.MEG = YGamma(n,:)';
    lme_gamma = fitlme(X,lme_formula); 
    gamma_coeff{n} = [lme_gamma.Coefficients.tStat(2); ...
        lme_gamma.Coefficients.pValue(2)];
    
    clc; fprintf('Done roi %.0f/%.0f\n',n,nrois);

end

aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

gamma_coeff = cell2mat(gamma_coeff);
alpha_coeff = cell2mat(alpha_coeff);
beta_coeff = cell2mat(beta_coeff);
theta_coeff = cell2mat(theta_coeff);
delta_coeff = cell2mat(delta_coeff);


%%

pv = 0.05;
M = nrois;
tlim = 6;

C = corr(YDelta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

figure; set(gcf,'color','w','name',lme_formula) 

subplot(2,3,1); 
semilogy(delta_coeff(2,:),'k'); 
hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(delta_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
xlabel('ROI');
title('Delta power'); grid on

aal_labels(delta_coeff(2,:)<=alpha)


C = corr(YTheta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);


subplot(2,3,2); 
semilogy(theta_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(theta_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Theta power'); grid on
xlabel('ROI');

aal_labels(theta_coeff(2,:)<=alpha)

C = corr(YAlpha');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,3); 
semilogy(alpha_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(alpha_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Alpha power'); grid on
xlabel('ROI');

aal_labels(alpha_coeff(2,:)<=alpha)

C = corr(YBeta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,4);
semilogy(beta_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r');
ylabel('p-value')
yyaxis right; plot(beta_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Beta power'); grid on
xlabel('ROI');

aal_labels(beta_coeff(2,:)<=alpha)

C = corr(YGamma');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,5); 
semilogy(gamma_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r');
ylabel('p-value')
yyaxis right; plot(gamma_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Gamma power'); grid on
xlabel('ROI');

aal_labels(gamma_coeff(2,:)<=alpha)