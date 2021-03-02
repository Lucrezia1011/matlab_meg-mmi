clear all
close all
% clc
meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

data_list = [];


% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
subexclude = [10,24];

roiopt = 'grid'; % running for grid
switch roiopt
    case 'grid'
        subexclude = [subexclude,26,49,53];
end

Nlist(subexclude) = []; 
zz= 0;
for sn = Nlist %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sdan = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;
            data_list{zz} = data_name;
        end
    end
    
   
end

%%
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca
% 
% % roiopt = 'c';
% roiopt = 'grid';
% % tfsopt = 'pwelch';
% tfsopt = 'm';
% 
gridres= 5;
% mu = 0.002; %saved as mu 0 
mu = 0.05; %mu=0.01
% freql=[ 4 8; 8 13; 13 25; 25 40; 40 150]; filter_type = 'but';
% freql = [1 4]; filter_type = 'firls';

freqband = [25 40]; 
% freqband = [40 100]; 
for ii = [1:length(data_list)]
    data_name = data_list{ii};
    sub = data_name(5:9);
%     mmiPreMoodPower(data_name,roiopt,gridres,freqband,mu); % mu = 0.05 
    mmiPreMoodPower_multiSpheres(data_name,roiopt,gridres,freqband,mu); % mu = 0.05 
end

%
% 
% if strcmp(roiopt,'grid') || strcmp(roiopt,'sens')
%     for freq = 1:size(freql,1)
%         for ii = 1:length(data_list)       
%             mmiPreMoodPower(data_list{ii},roiopt,gridres,freql(freq,:)); % pre mood
% %             mmi_grid_prep_PowerITI(param_list{ii},roiopt,gridres,freql(freq,:))
% %             mmi_grid_prep_Power(param_list{ii},roiopt,gridres,freql(freq,:)) % pre mood
% %             mmi_grid_prep_PowerA(param_list{ii},roiopt,gridres,freql(freq,:),filter_type)
% %             % anticipation period          
%         end
%     end    
% else
%     for ii = 1:length(data_list)
%         mmi_grid_prep(data_list{ii},roiopt,gridres,tfsopt)
%     end
% end

% resave ltv variables from delta

return

%%
Pall = [];

ltvMood = [];
% 11,12 (sub 22695)   58,59 sub(23927)
for sn = 1:length(data_list)
    
    
    data_name = data_list{sn};
    sub = data_name(5:9);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
    
    
    
    cd(data_path)
    clear grid
    load('leadfields_multiSpheres_5mm.mat');
    if nnz(grid.inside) == 0
        data_name
    end
    if sn == 1
        gridall = zeros(length(grid.inside),length(data_list));
    end
     gridall(:,sn)  = grid.inside;
    
    
    load(sprintf('pre_mood_grid_%.0f-%.0fHz_mu%g_multiSpheres',...
        freqband(1),freqband(2),mu*100))
    %         load(sprintf('pre_mood_grid_%.0f-%.0fHz',...
    %             freqband(1),freqband(2)))
%     P = (P - mean(P(:))) / std(P(:));
    Pmni = zeros(size(grid.inside,1),size(P,2));
    Pmni(grid.inside,:) = P;
    Pall = cat(2,Pall,Pmni);
    
    
    ltvmood.recording = repmat(sn,size(ltvmood,1),1);
    if isempty(ltvMood)
        ltvMood = ltvmood;
    else
        ltvMood(end+(1:size(ltvmood,1)),:) = ltvmood;
    end
    
    
end

gridall = all(gridall,2);
clear ltvmood TFS P Pmni grid
if strcmp(roiopt,'grid')
    Pall = Pall(gridall,:); 
end


%% Write data (grid)
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';
% dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
%     freqband(1),freqband(2)),Pall);
dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz_mu%g_multiSpheres.txt',out_path,...
    freqband(1),freqband(2),mu*100),Pall);
dlmwrite([out_path,'/mni_grid_multiSpheres.txt'],gridall);
% writetable(ltvMood,[out_path,'/latent_vars.csv']);
 

%%

% mmiPreMoodPower_plotGrid.m
% Lucrezia Liuzzi, last edited 2020/12/02
% Plot Linear mixed effects model results for high-beta power on MNI brain

% start fieldtrip
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';
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
gridall = dlmread([out_path,'/mni_grid_multiSpheres.txt']);


%% Plot significant results

datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz_mu5_multiSpheres',...
    freql(1),freql(2));
muname = 'mu5';

nf = size(freql,1);
Tmood = zeros(size(gridall));
Te = zeros(size(gridall));
Tesum = zeros(size(gridall));

param_list = cell(1,21);
for nn = 1:21
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end


cd([datapath,'/lme_E_sum'])

X = zeros(nnz(gridall),1);
for nn = 1:21
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
end

Tesum(gridall==1) = X;


cd([datapath,'/lme_E'])

X = zeros(nnz(gridall),1);
pV =  zeros(nnz(gridall),1);
for nn = 1:21
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
end
Te(gridall==1) = X;

fprintf(sprintf('Loaded frequency %.0f-%.0fHz\n',freql(1),freql(2)))

%% Read random permutations
% Declared in pre-registratiob to use 2-tailed t-test with 2,000 randperms
% Ranperms for E and E_sum, mu = 0.05 only
 
fit_parameters =  {'E';'E_sum'};
   
% Read random permutations: t-statistic from fixed effect 
clusternull = cell(1,2);
for fitpar = 1 
    
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
%% Plot significant clusters

for p =  2
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

%     ft_sourceplot(cfg, sourceout_Int);
%     saveas(gcf,sprintf('~/matlab/figures/%s_betagamma_peak3_%s_multiSpheres_1e3perms.png',fit_parameter,muname))

    % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz\npeak t-value %.1f',...
    %     H,E,fit_parameter,freql(1),freql(2),max(abs(sourceant.pow(:)))))
    % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz',...
    %     H,E,fit_parameter,freql(1),freql(2)))
    cfg.method        = 'slice'; %'ortho'
    ft_sourceplot(cfg, sourceout_Int);
    saveas(gcf,sprintf('~/matlab/figures/%s_betagamma_slice_%s_multiSpheres_1e3perms.png',fit_parameter,muname))
end




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
