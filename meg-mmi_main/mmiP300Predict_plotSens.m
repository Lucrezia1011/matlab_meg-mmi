% Lucrezia Liuzzi, last updated 2021/03/25
% 30Hz highpassed MEG signal 250-400ms after gambling options presentation. 
% Plot results from sensor level analysis

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% Load data   
param_list{1} = '001';
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/P300/confirm/';

if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    sensall = channels;
end

freql = {'cue'};
ff = 1;

latent_vars_name = sprintf('latent_vars_%s.csv',freql{ff});
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
fit_parameters = X.Properties.VariableNames(5:6);

freq = sprintf('M%s_P300_30Hzlowpass',freql{ff});
% It does not work with planar gradiometers. Results completely different
% from model with original axial data. Planar gradiometer transform changes
% distribution of noise and cannot be used for trial-level analysis.
% freq = sprintf('M%s_P300_30Hzlowpass_planar',freql{ff});
% fit_parameters = X.Properties.VariableNames([5,6,8]);


meg_data_name = sprintf('%s.txt',freq);
meg = dlmread([data_path,meg_data_name]);

outpath = sprintf('%s%s/',data_path,freq);
nn =1;

Tfit = cell(1,length(fit_parameters));
pfit = cell(1,length(fit_parameters));
for ii = 1:length(fit_parameters)
    cd([outpath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    Tfit{ii} = Xv.tStat;
    pfit{ii} = Xv.pValue;
end

for ii = 1:2
    [mxt,mxii] = max(Tfit{ii});
    fprintf('Peak t-values for %s:\nt-value=%.2f, p-value=%e\n',fit_parameters{ii},mxt,pfit{ii}(mxii))
    [mxt,mxii] = min(Tfit{ii});
    fprintf('Peak t-values for %s:\nt-value=%.2f, p-value=%e\n',fit_parameters{ii},mxt,pfit{ii}(mxii))
end


T = struct;
T.label = sensall;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position',[176 348 1262 800])

for ii = 1:length(fit_parameters)
    subplot(2,3,ii)
    T.trial{1} = Tfit{ii}; T.avg = Tfit{ii};
    cfg = [];
    cfg.channel = sensall;
    cfg.layout = 'CTF275_helmet.mat';
    cfg.zlim = [-5 5];
    cfg.comment ='no';
    ft_topoplotER(cfg, T)
    colorbar
    titlename = fit_parameters{ii};
    k = strfind(titlename,'_');
    titlename(k) = ' ';
    title(titlename)
end
% saveas(gcf,sprintf('~/matlab/figures/%s.png',freq))


% plot the MEG data average(z-score)
subs = unique(X.subject);
megz = zeros(size(meg,1),length(subs));
for ss = 1:length(subs)
    megm = mean(meg(:,X.subject==subs(ss)),2); % mean over trials
    megz(:,ss)  = zscore(megm);
end
figure; clf; set(gcf,'color','w')
T.trial{1} = mean(megz,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-.8 .8 ];
cfg.colorbar = 'yes';
cfg.comment = 'no';
ft_topoplotER(cfg, T);
title(sprintf('Grand-average (Z-scored)'))

% plot the MEG data average
figure; set(gcf,'color','w')
T.trial{1} = mean(meg,2); T.avg = mean(meg,2);
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-1 1]*3e-14;
cfg.comment ='no';
ft_topoplotER(cfg, T)
titlename = 'Grandaverage of MEG sensor data';
colorbar
% saveas(gcf,sprintf('~/matlab/figures/%s_grandav.png',freq))

%% Read randperms
clusternull = cell(1,2);

cd([outpath,'/lme_E_LTA/'])
% 
% nullnames = dir('sens_permute');
% nullnames(1:2)=[];
% % Delete empty data
% for n = 1:length(nullnames)
%     if nullnames(n).bytes == 0
%         delete(['sens_permute/',nullnames(n).name])
%     end
% end
nullnames = dir('sens_permute');
nullnames(1:2)=[];

clusternull{1} = zeros(length(nullnames),size(meg,1));
for n = 1:length(nullnames)
    clusternull2 = dlmread(['sens_permute/',nullnames(n).name]);
    clusternull{1}(n,:) = clusternull2;
end
figure; histogram(clusternull{1}(:))


cd([outpath,'/lme_E_sum/'])

% nullnames = dir('sens_permute');
% nullnames(1:2)=[];
% % Delete empty data
% for n = 1:length(nullnames)
%     if nullnames(n).bytes == 0
%         delete(['sens_permute/',nullnames(n).name])
%     end
% end

nullnames = dir('sens_permute');
nullnames(1:2)=[];

clusternull{2} = zeros(length(nullnames),size(meg,1));
for n = 1:length(nullnames)
    clusternull2 = dlmread(['sens_permute/',nullnames(n).name]);
    clusternull{2}(n,:) = clusternull2;
end
hold on; histogram(clusternull{2}(:))

%% Plot p-values with multiple comparison correction

figure; set(gcf,'color','w','position', [122   460   828   387])

for ii = 1:2
    subplot(1,2,ii)
    p = sort(pfit{ii});
    semilogy(p)
    hold on
    semilogy(0.05./(length(sensall):-1:1))
%     grid on
    xlabel('sensors')
    ylabel('p-value')
    legend('p-values','FDR','location','best')
    N = nnz(p'<0.05./(length(sensall):-1:1));
    
    titlename = fit_parameters{ii};
    k = strfind(titlename,'_');
    titlename(k) = ' ';
    title(sprintf('%s: p-value of %.0f sensors < 0.05 (FDR)',titlename,N))
end
% saveas(gcf,sprintf('~/matlab/figures/%s_pvalue.png',freq))

%% Clustering of random permutations
% Read data header from one subject to obtain sensors positions
cd(['/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/'])
hdr = ft_read_header('sub-24071_task-mmi3_run-1_meg.ds');

cfg_neighb        = [];
cfg_neighb.method = 'template';%'distance';
cfg_neighb.channel = channels;
cfg_neighb.template = 'CTF275_neighb.mat';
neighbours        = ft_prepare_neighbours(cfg_neighb, hdr);
for n = 1:length(neighbours)
    [~,~,iB] = intersect(neighbours(n).neighblabel, channels );
    neighbours(n).neighbnum =iB;
end


% Declared in pre-reg to use 2-tailed t-test with 2,000 randperms
E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1;
% Permutations of the linear mixed model fixed effect
clusternull = cell2mat(clusternull');
% M = zeros(2,length(clusternull));
M = cell(1,length(clusternull));

parfor n = 1:length(clusternull)
    
    img = clusternull(n,:);    
    [tfced] = matlab_tfce_transform_MEG(img,H,E,dh,neighbours);
    [tfcen] = -matlab_tfce_transform_MEG(-img,H,E,dh,neighbours);

%     M(1,n) = min(tfcen);
%     M(2,n) = max(tfced);
    M{n} = [min(tfcen), max(tfced)];
    if mod(n,200) == 0
        clc; fprintf('TFCE null done %.0f perc.\n',n/length(clusternull)*100)
    end
end
%% Plots
% M = cell2mat(M);
alpha = 0.05; 
alpha = alpha/2; % two-tailed
nparam = 2; % number of predictors for multiple comparisons
Mp = sort(M,'descend');
threshp = Mp(round(alpha/nparam*size(clusternull,1)));
Mn = sort(M,'ascend');
threshn = Mn(round(alpha/nparam*size(clusternull,1)));

figure; clf; set(gcf,'color','w','position',[176 748 1262 400])
c0 = get(gca,'colormap');
for ii = 1:2
[tfced] = matlab_tfce_transform_MEG(Tfit{ii}',H,E,dh,neighbours);
[tfcen] = matlab_tfce_transform_MEG(-Tfit{ii}',H,E,dh,neighbours);
tfced = tfced - tfcen;
Tplot = Tfit{ii};
Tplot(tfced<threshp & tfced>threshn) = 0;

clim = [-5,5];
a = linspace(clim(1),clim(2),256);
c = c0;
cthresh = min(Tplot(Tplot>0));
c(abs(a)<2.5,:) = repmat([1 1 1],nnz(abs(a)<2.5),1);

T = struct;
T.label = channels;
T.time{1} = 300;
T.sampleinfo = [1 1];
T.trial{1} =Tplot ;
T.avg = T.trial{1};% T.mask = sigS(sind)';

cfg = [];
% cfg.mask = Tplot; % masking doesn't work
% cfg.maskparameter = 'mask';
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
cfg.colormap   = c;
cfg.zlim =clim;
cfg.comment    = 'no';
subplot(1,2,ii)
ft_topoplotER(cfg, T)
colorbar
end

saveas(gcf,sprintf('~/matlab/figures/P300_sens_1e4perms_confirm.png'))

%% Get peak sensor as well as whole cluster
for ii = 1:2
[tfced] = matlab_tfce_transform_MEG(Tfit{ii}',H,E,dh,neighbours);
[tfcen] = matlab_tfce_transform_MEG(-Tfit{ii}',H,E,dh,neighbours);
tfced = tfced - tfcen;
Tplot = Tfit{ii};
Tplot(tfced<threshp & tfced>threshn) = 0;
C = channels(Tplot~=0);
save(sprintf('%s/sig_channels_%s.mat',data_path,fit_parameters{ii}),'C')
end
