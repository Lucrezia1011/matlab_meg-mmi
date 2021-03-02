% mmiPreMoodPower_plotGrid.m
% Lucrezia Liuzzi, last edited 2020/12/02
% Plot Linear mixed effects model results for high-beta power on MNI brain

% start fieldtrip
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/confirm/';
ltvMood = readtable([out_path,'/latent_vars.csv']);

freql=[ 25 40]; % frequency band

% Parameters for TFCE
E = 0.5;
H = 2; 
dh = 0.1;

% Load standard brain for plotting
mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
% mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates


ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5; % grid resolution 5mm
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

% MNI voxels included in model
gridall = dlmread([out_path,'/mni_grid.txt']);

%% Read random permutations
% Declared in pre-registratiob to use 2-tailed t-test with 2,000 randperms
% Ranperms for E and E_sum, mu = 0.05 only

datapath = sprintf('%s/powergrid_%.0f-%.0fHz_mu5',...
            out_path,freql(1),freql(2));
 
fit_parameters =  {'E';'E_sum'};
   
% Read random permutations: t-statistic from fixed effect 
clusternull = cell(1,2);
for fitpar = 1:2  
    
    cd([datapath,'/lme_',fit_parameters{fitpar},'/'])

    nullnames = dir('grid_permute');
    nullnames(1:2)=[];
    % Delete empty data
    for n = 1:length(nullnames)
        if nullnames(n).bytes == 0
            delete(['grid_permute/',nullnames(n).name])
        end
    end

    nullnames = dir('grid_permute');
    nullnames(1:2)=[];

    clusternull{fitpar} = zeros(length(nullnames),nnz(gridall));
    for n = 1:length(nullnames)
        clusternull2 = dlmread(['grid_permute/',nullnames(n).name]);
        clusternull{fitpar}(n,:) = clusternull2;
    end
end

figure(1); clf; subplot(121);
histogram(clusternull{1}(:),'DisplayStyle','stairs')
hold on; histogram(clusternull{2}(:),'DisplayStyle','stairs')
title('t-stat from random permutations (all voxels)')
xlabel('t-stat'); ylabel('no. voxels (permutations x mni voxels)')
legend(fit_parameters)

clusternull = cell2mat(clusternull');
M = cell(1,length(clusternull));

% Create null ditribution from max/min TFCE over all brain
parfor n = 1:length(nullnames)
    clusternull2 = zeros(size(gridall));
    clusternull2(gridall==1) = clusternull(n,:);
    img = reshape(clusternull2,sourcemodel.dim);
    
    tfcep= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
    
    tfcen= -matlab_tfce_transform(-img,H,E,26,dh); % C=26 default setting
    
    M{n} = [min(tfcen(:)) max(tfcep(:))];
end


M = cell2mat(M);
alpha = 0.05; 
alpha = alpha/2; % two-tailed
nparam = 2; % number of predictors for multiple comparisons 'E' and 'Esum'
Mp = sort(M,'descend');
threshp = Mp(round(alpha/nparam*size(clusternull,1)));
Mn = sort(M,'ascend');
threshn = Mn(round(alpha/nparam*size(clusternull,1)));

% Plot null distribution of TFCE values
figure(1); subplot(122); 
histogram(M)
hold on
y = ylim;
plot(threshn*[1 1],[0 y(2)],'r')
plot(threshp*[1 1],[0 y(2)],'r')
title('t-stat from random permutations (all voxels)')
xlabel('TFCE'); ylabel('no. permutations')
ylabel('TFCE null distribution')
%% Plot significant results

% loop over 3 regularization parameters
for mu = 4
switch mu
    case 1
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz',...
            freql(1),freql(2));
        muname = 'mu5';
    case 2
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu1',...
            freql(1),freql(2));
        muname = 'mu1';
    case 3
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu0.2',...
            freql(1),freql(2));
        muname = 'mu0.2';
    case 4 % confirmatory analysis
        datapath = sprintf('%s/powergrid_%.0f-%.0fHz_mu5',...
            out_path,freql(1),freql(2));
        muname = 'mu5';
end

nf = size(freql,1);
Tmood = zeros(size(gridall));
Te = zeros(size(gridall));
Tesum = zeros(size(gridall));
Pe = zeros(size(gridall));
Pesum = zeros(size(gridall));

param_list = cell(1,15);
for nn = 1:15
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

% cd([datapath,'/lme_mood'])
% cd([datapath,'/lme_M'])
% 
% X = zeros(nnz(gridall),1);
% for nn = 1:14
%     opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
%     Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
%     X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
% end
% Tmood(gridall==1) = X;
% 
% 
cd([datapath,'/lme_E_sum'])

X = zeros(nnz(gridall),1);
pV =  zeros(nnz(gridall),1);

for nn = 1:15
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
end

fprintf('Fixed effect E_sum:\nPeak T-stat=%.2f,  p-value=%.1e\n',...
        max(X),min(pV))
Tesum(gridall==1) = X;
Pesum(gridall==1) = pV;

cd([datapath,'/lme_E'])

X = zeros(nnz(gridall),1);
pV =  zeros(nnz(gridall),1);
for nn = 1:15
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
end
fprintf('Fixed effect E_LTA:\nPeak T-stat=%.2f,  p-value=%.1e\n',...
        max(X),min(pV))
Te(gridall==1) = X;
Pe(gridall==1) = pV;

fprintf(sprintf('Loaded frequency %.0f-%.0fHz\n',freql(1),freql(2)))


%% Plot significant clusters

for p =  2:3
    sourceant =[];

    switch p
        case 1
            sourceant.pow = Tmood;
            fit_parameter = 'mood';
        case 2
            sourceant.pow = Te;
            fit_parameter = 'E';
        case 3
            sourceant.pow = Tesum;
            fit_parameter = 'E sum';
    end
      
    % Use spatial permutations to find clusters
    % Smith and Nichols, Threshold-Free Cluster Enhancement, 2009. 
    % "For example, to correct for multiple comparisons, one simply has to build
    % up the null distribution (across permutations of the input data) of the 
    % maximum (across voxels) TFCE score, and then test the actual TFCE image 
    % against that. Once the 95th percentile in the null distribution is found 
    % then the TFCE image is simply thresholded at this level to give inference 
    % at the p<0.05 (corrected) level."
    
%     nperms = 1e4;
%     M = zeros(1,length(nullnames));
%     T = sourceant.pow(gridall==1);
%     for n = 1:nperms
%         clusternull2 = zeros(size(gridall));
%         clusternull2(gridall==1) = T(randperm(nnz(gridall)));
% 
%         img = reshape(clusternull2,sourcemodel.dim);
% 
%         tfce = matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
%         tfce = sort(tfce(:),'descend');
%         M(n) = max(tfce(:));
%         if mod(n,100)==0
%             clc
%             fprintf('Done %.0f perc.\n',n/nperms*1e2)
%         end
%     end
%     M = sort(M,'descend');
%     thresh = M(round(0.05*nperms)); %mu=0.2% thresh=796.9783
    
    
    img = reshape(sourceant.pow,sourcemodel.dim);

    tfce= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
    tfcen= matlab_tfce_transform(-img,H,E,26,dh); % C=26 default setting
    tfce = tfce - tfcen;
    sourceant.pow(tfce>threshn & tfce<threshp)  = 0;
    % sourceant.pow  = tfce(:);
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;

    imgtfce = reshape(sourceant.pow,sourceant.dim);
    peaks3d = imregionalmax(imgtfce);
    [~,peaks3dSort] = sort(sourceant.pow(peaks3d),'descend');
    peaks3dind = find(peaks3d);
    peaks3dcoord = sourceant.pos(peaks3dind(peaks3dSort),:)*10;

    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'mni';

    imgtfce = sourceout_Int.pow;
    peaks3d = imregionalmax(imgtfce);
    [~,peaks3dSort] = sort(sourceout_Int.pow(peaks3d),'descend');
    peaks3dind = find(peaks3d);
    % aal.tissuelabel{aal.tissue(peaks3dind(peaks3dSort))}

    crang = [2,5];
    % crang = [thresh max(sourceant.pow)];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
    if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
        cfg.location   = 'max';
    else
        cfg.location   = 'min';
    end
%     cfg.location = [-2 -40 30]; % peak for mu =0.05
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

    ft_sourceplot(cfg, sourceout_Int);
    saveas(gcf,sprintf('~/matlab/figures/%s_gamma_peak_%s_2e3perms_confirm.png',fit_parameter,muname))

    % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz\npeak t-value %.1f',...
    %     H,E,fit_parameter,freql(1),freql(2),max(abs(sourceant.pow(:)))))
    % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz',...
    %     H,E,fit_parameter,freql(1),freql(2)))
    cfg.method        = 'slice'; %'ortho'
    ft_sourceplot(cfg, sourceout_Int);
    saveas(gcf,sprintf('~/matlab/figures/%s_gamma_slice_%s_2e3perms_confirm..png',fit_parameter,muname))
end
% close all
end


%% Plot source data from 1 recording
% megP = dlmread([datapath,'.txt']);
% 
% for ii = 1:max(ltvMood.recording)
%     M = megP(:,ltvMood.recording==ii);
%     megP(:,ltvMood.recording==ii) = (M) / std(M(:));
% end
% meg = zeros(size(sourcemodel.inside));
% meg(gridall==1) = mean(megP,2); % simple mean
% 
% ii = 6;
% 
% meg = zeros(size(sourcemodel.inside));
% M = megP(:,ltvMood.recording==ii);
% meg(gridall==1) = mean(M,2); 
% 
% sourceant.pow = meg;
% sourceant.dim = sourcemodel.dim;
% sourceant.inside = sourcemodel.inside;
% sourceant.pos = sourcemodel.pos;
% cfg = [];
% cfg.parameter = 'pow';
% sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
% sourceout_Int.pow(~sourceout_Int.inside) = 0;
% sourceout_Int.coordsys = 'mni';
% 
% 
% crang = [];
% cfg = [];
% cfg.method        = 'slice'; %'ortho'
% if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
%     cfg.location   = 'max';
% else
%     cfg.location   = 'min';
% end
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
% 
% ft_sourceplot(cfg, sourceout_Int);


%% Plot average power
for mu = [1,4]%1:4
switch mu
    case 1
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz',...
            freql(1),freql(2));
        muname = 'mu5';
        gridall = dlmread([out_path,'/mni_grid.txt']);
    case 2
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu1',...
            freql(1),freql(2));
        muname = 'mu1';
    case 3
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu0.2',...
            freql(1),freql(2));
        muname = 'mu0.2';
    case 4
        datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu5_multiSpheres',...
            freql(1),freql(2));
        muname = 'mu5';
        gridall = dlmread([out_path,'/mni_grid_multiSpheres.txt']);
end

    P = dlmread([datapath,'.txt']);
    subs = unique(ltvMood.subject);
    Ps = zeros(size(P,1),length(subs));
    for s = 1:length(subs)
        Ps(:,s) = (mean(P(:,ltvMood.subject == subs(s,:)),2));
    end

    sourceant =[];
    sourceant.pow = zeros(size(gridall));
    sourceant.pow(gridall==1) = mean(Ps,2);
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    sourceout_Int.pow(~sourceout_Int.inside) = 0;
    sourceout_Int.coordsys = 'mni';

    crang = [];
%     crang = [0,max(sourceout_Int.pow(:))];
    cfg = [];
    cfg.method        = 'ortho'; %'ortho'
    if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
        cfg.location   = 'max';
    else
        cfg.location   = 'min';
    end
%     cfg.location = [-2 -40 30]; % peak for mu =0.05
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

    ft_sourceplot(cfg, sourceout_Int);
    title(muname)
%     saveas(gcf,sprintf('~/matlab/figures/gamma_av_%s.png',muname))
end