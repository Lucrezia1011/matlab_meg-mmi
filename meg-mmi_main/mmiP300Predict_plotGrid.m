% Lucrezia Liuzzi, last updated 2021/03/25
% 30Hz highpassed MEG signal 250-400ms after gambling options presentation. 
% Plot results from source level analysis

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% Load data and set parameters

mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/P300/confirm/';
freql = {'cue';'choice'}; 
ff = 1;
latent_vars_name = sprintf('latent_vars_%s.csv',freql{ff});
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
fit_parameters = X.Properties.VariableNames(5:6);

% gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/mni_grid.txt');
gridall = dlmread([data_path,'mni_grid.txt']);

% Parameters for TFCE
E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1;

alpha = 0.05;


%% Read random permutations
% Declared to use 2-tailed t-test with 2,000 randperms
% Ranperms for E and E_sum, mu = 0.05 only

datapath = sprintf('%sBF%s_P300_30Hzlowpass/',data_path,freql{ff});
 
clusternull = cell(1,2);
for fitpar = 1:2
    
    cd([datapath,'/lme_',fit_parameters{fitpar},'/'])

%     nullnames = dir('grid_permute');
%     nullnames(1:2)=[];
    % Delete empty data
%     for n = 1:length(nullnames)
%         if nullnames(n).bytes == 0
%             delete(['grid_permute/',nullnames(n).name])
%         end
%     end

    nullnames = dir('grid_permute');
    nullnames(1:2)=[];

    clusternull{fitpar} = zeros(length(nullnames),nnz(gridall));
    for n = 1:length(nullnames)
        clusternull2 = dlmread(['grid_permute/',nullnames(n).name]);
        clusternull{fitpar}(n,:) = clusternull2;
    end
end

figure(1); clf; subplot(121);
histogram(clusternull{1},'DisplayStyle','stairs','Normalization','pdf')
hold on; histogram(clusternull{2}(:))

% Permutations of the linear mixed model fixed effect
clusternull = cell2mat(clusternull');
M = cell(1,size(clusternull,1));

parfor n = 1:size(clusternull,1)
    clusternull2 = zeros(size(gridall));
    clusternull2(gridall==1) = clusternull(n,:);
    img = reshape(clusternull2,sourcemodel.dim);
    
    tfcep= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
    
    tfcen= -matlab_tfce_transform(-img,H,E,26,dh); % C=26 default setting
    
    M{n} = [min(tfcen(:)) max(tfcep(:))];
end


M = cell2mat(M);

alpha = alpha/2; % two-tailed
nparam = 2; % number of predictors for multiple comparisons
Mp = sort(M,'descend');
threshp = Mp(round(alpha/nparam*size(clusternull,1)));
Mn = sort(M,'ascend');
threshn = Mn(round(alpha/nparam*size(clusternull,1)));

figure(1); subplot(122); 
histogram(M)
hold on
plot(threshn*[1 1],[0 size(clusternull,1)/2],'r')
plot(threshp*[1 1],[0 size(clusternull,1)/2],'r')


%% Read results of mixed effects model
clear T pV
T(1:length(fit_parameters)) = {zeros(size(gridall))};
pV(1:length(fit_parameters)) = {zeros(size(gridall))};
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
    
for ii = 1:length(fit_parameters)
    cd([datapath,'lme_',fit_parameters{ii}])
    X = zeros(nnz(gridall),1);
    pv = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
        pv((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
    end
    T{ii}(gridall==1) = X;
    pV{ii}(gridall==1) = pv;
    
end


for ii = 1:2
    [mxt,mxii] = min(T{ii}); % peaks are negative
    fprintf('Peak t-values for %s:\nt-value=%.2f, p-value=%e\n',fit_parameters{ii},mxt,pV{ii}(mxii))
end


%% Plot pValues with FDR correction
freq = sprintf('BF%s_P300_30Hzlowpass',freql{ff});
% meg_data_name = sprintf('%s.txt',freq);
% meg = dlmread([data_path,meg_data_name]);
% 
% M = size(meg,1);
% 
% C = corr(meg');
% lambda = eig(C);
% % Effective number of independent variables
% Meff = 1 + (M-1)*(1 - var(lambda)/M);
% alpha = 1 - (1 - 0.05)^(1/Meff);


figure; set(gcf,'color','w')

for ii = 1:length(fit_parameters)
    
    subplot(2,3,ii)
    
    p = sort(pV{ii}(gridall==1));
    semilogy(p)
    hold on
    semilogy(0.05./(nnz(gridall):-1:1))
    xlim([1 100])
%     grid on
    xlabel('voxels')
    ylabel('p-value')
    legend('p-values','FDR','location','best')
    N = nnz(p'<0.05./(nnz(gridall):-1:1));
    
    title(sprintf('%s: p-value of %.0f voxels < 0.05 (FDR)',fit_parameters{ii},N))
end


%% Plot significant clusters

for ii = 1:2
    sourceant =[];
    sourceant.pow = T{ii};
    fit_parameter = fit_parameters{ii};
    
      
    % Use spatial permutations to find clusters
    % Smith and Nichols, Threshold-Free Cluster Enhancement, 2009. 
    % "For example, to correct for multiple comparisons, one simply has to build
    % up the null distribution (across permutations of the input data) of the 
    % maximum (across voxels) TFCE score, and then test the actual TFCE image 
    % against that. Once the 95th percentile in the null distribution is found 
    % then the TFCE image is simply thresholded at this level to give inference 
    % at the p<0.05 (corrected) level."
 
    img = reshape(sourceant.pow,sourcemodel.dim);

    tfce= matlab_tfce_transform(img,H,E,26,dh); % C=26 default setting
    tfcen= matlab_tfce_transform(-img,H,E,26,dh); % C=26 default setting
    tfce = tfce - tfcen;
    sourceant.pow(tfce>threshn & tfce<threshp)  = 0;
    
    if all(sourceant.pow==0)
        fprintf('%s is not significant\n',fit_parameter)
    else
        
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

        crang = [];
        % crang = [thresh max(sourceant.pow)];
        cfg = [];
        cfg.method        = 'ortho'; %'ortho'
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
        saveas(gcf,sprintf('~/matlab/figures/%s_P300_peak_2e3perms_confirm.png',fit_parameter))

        % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz\npeak t-value %.1f',...
        %     H,E,fit_parameter,freql(1),freql(2),max(abs(sourceant.pow(:)))))
        % title(sprintf('H=%.1f, E=%.1f .  %s : %.0f-%.0fHz',...
        %     H,E,fit_parameter,freql(1),freql(2)))
        cfg.method        = 'slice'; %'ortho'
        ft_sourceplot(cfg, sourceout_Int);
        saveas(gcf,sprintf('~/matlab/figures/%s_P300_slice_2e3perms_confirm.png',fit_parameter))
    end
end

%% Plot timecourse


data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/';
Y = dlmread([data_path,'/meg_trials_evoked_cue.txt']);
Y = reshape(Y, 360,116,size(Y,2));
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');
rois = 104; %peak on Cerebellum 8 R

Nlist = 1:size(Y,3);
% for confirmatory analysis
subexclude =[1:12,14:16]; %but 10 is already excluded from Y
subexclude(10) = [];
subexclude(subexclude>10) =subexclude(subexclude>10)-1;
Nlist(subexclude) = []; 
Y = Y(:,:,Nlist);

time = linspace(-.2,1,360);
figure(3);  clf;set(gcf,'color','w')
plot(time,mean(Y(:,rois,:),3))
hold on
fill([0.25 0.4 .4 0.25],[-1 -1 1 1]*.4,[0.5 0.5 0.5],'facealpha',0.3,'edgecolor','none')

fill([time fliplr(time)],[mean(Y(:,rois(1),:),3)-std(Y(:,rois(1),:),0,3)/sqrt(nnz(Nlist)); ...
    flipud(mean(Y(:,rois(1),:),3)+std(Y(:,rois(1),:),0,3)/sqrt(nnz(Nlist)))],[0 0 1],'facealpha',0.3,'edgecolor','none')

% fill([time fliplr(time)],[mean(Y(:,rois(2),:),3)-std(Y(:,rois(2),:),0,3)/sqrt(nnz(Nlist)); ...
%     flipud(mean(Y(:,rois(2),:),3)+std(Y(:,rois(2),:),0,3)/sqrt(nnz(Nlist)))],[1 0 0],'facealpha',0.3,'edgecolor','none')

grid on
% ylim([-1 1]*0.6e2);
xlabel('time (s)'); ylabel('Beamformed evoked response (z-score)')
title(sprintf('Average signal (no. subjects=%d) from significant ROIs',nnz(Nlist)))
legend(aal_labels{rois(1)},...
    'Hypothesis time window',...
    'Standard error over subjects',...
    'location','best') % standard error over subjects

saveas(gcf,sprintf('~/matlab/figures/P300_sourcepeak_timecourse.png'))