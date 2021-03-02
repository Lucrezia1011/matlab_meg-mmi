%%  Plot Linear mixed effects model for sensors
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

param_list{1} = '001';
freql=[1 4; 4 8; 8 13; 13 25; 25 40; 50 150];
if ~exist('sensall','var')
%     load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens//pre_mood/confirm/sensors.mat')
end

% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/';
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/confirm/';

ff = 5; %ff = 5;


freq = sprintf('powersens_%.0f-%.0fHz',...
    freql(ff,1),freql(ff,2));

meg_data_name = sprintf('%s.txt',freq);
meg = dlmread([data_path,meg_data_name]);

outpath = sprintf('%s%s/',data_path,freq);
nn =1;
% fit_parameters = {'mood';'E';'E_sum';'RPE';'RPE_sum'};
fit_parameters = {'E';'E_sum'};
Tfit = cell(1,length(fit_parameters));
pfit = cell(1,length(fit_parameters));
for ii = 1:length(fit_parameters)
    cd([outpath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    Tfit{ii} = Xv.tStat;
    pfit{ii} = Xv.pValue;
    
    [~,p] = min(Xv.pValue);
    fprintf('Fixed effect %s:\nPeak on channel %s\nT-stat=%.2f,  p-value=%.1e\n',...
        fit_parameters{ii},channels{p},Xv.tStat(p),Xv.pValue(p))
end


%% Read randperms
clusternull = cell(1,2);

cd([outpath,'/lme_E/'])

nullnames = dir('sens_permute');
nullnames(1:2)=[];
% Delete empty data
for n = 1:length(nullnames)
    if nullnames(n).bytes == 0
        delete(['sens_permute/',nullnames(n).name])
    end
end
nullnames = dir('sens_permute');
nullnames(1:2)=[];

clusternull{1} = zeros(length(nullnames),size(meg,1));
for n = 1:length(nullnames)
    clusternull2 = dlmread(['sens_permute/',nullnames(n).name]);
    clusternull{1}(n,:) = clusternull2;
end
figure; histogram(clusternull{1}(:))


cd([outpath,'/lme_E_sum/'])

nullnames = dir('sens_permute');
nullnames(1:2)=[];
% Delete empty data
for n = 1:length(nullnames)
    if nullnames(n).bytes == 0
        delete(['sens_permute/',nullnames(n).name])
    end
end

nullnames = dir('sens_permute');
nullnames(1:2)=[];

clusternull{2} = zeros(length(nullnames),size(meg,1));
for n = 1:length(nullnames)
    clusternull2 = dlmread(['sens_permute/',nullnames(n).name]);
    clusternull{2}(n,:) = clusternull2;
end
hold on; histogram(clusternull{2}(:))

%%

T = struct;
T.label = channels;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position', [387  558  1139  416])

for ii = 1:length(fit_parameters)
subplot(1,3,ii)
T.trial{1} = Tfit{ii}; T.avg = Tfit{ii};
cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-5 5];
ft_topoplotER(cfg, T)

titlename = fit_parameters{ii};
k = strfind(titlename,'_');
titlename(k) = ' ';
title(titlename)

end

% saveas(gcf,sprintf('~/matlab/figures/pre_mood_%s.png',freq))

% plot the MEG data
ltvMood = readtable([data_path,'/latent_vars.csv']);
subs = unique(ltvMood.subject);
megz = zeros(size(meg,1),length(subs));
for ss = 1:length(subs)
    megm = mean(meg(:,ltvMood.subject==subs(ss)),2); % mean over trials
    megz(:,ss)  = zscore(megm);
end
figure(ff+2); clf; set(gcf,'color','w','position', [387  558  1139  416])
subplot(121)
T.trial{1} = mean(meg,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [0.5 2]*1e-27;
cfg.colorbar = 'yes';
cfg.comment = 'no';
ft_topoplotER(cfg, T);
title(sprintf('Grand-average of 25-40Hz variance\nduring mood rating preparation'))

subplot(122)
T.trial{1} = mean(megz,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-1 1];
cfg.colorbar = 'yes';
cfg.comment = 'no';
ft_topoplotER(cfg, T);
title(sprintf('Grand-average (Z-scored) of 25-40Hz variance\nduring mood rating preparation'))
saveas(gcf,sprintf('~/matlab/figures/pre_mood_%s_avg.png',freq))

% Plot p-values with multiple comparison correction

figure; set(gcf,'color','w')

ii = 2;

p = sort(pfit{ii});
semilogy(p)
hold on
semilogy(0.05./(length(channels):-1:1))
%     grid on
xlabel('sensors')
ylabel('p-value')
legend('p-values','FDR','location','best')
N = nnz(p'<0.05./(length(channels):-1:1));

title(sprintf('%s: p-value of %.0f sensors < 0.05 (FDR)',fit_parameters{ii},N))

% saveas(gcf,sprintf('~/matlab/figures/pre_mood_%s_pvalue.png',freq))

% 
% T = struct;
% T.label = channels;
% T.time{1} = 300;
% T.sampleinfo = [1 1];
% figure; clf; set(gcf,'color','w','position',[176 748 1262 400])
% 
% c = get(gca,'colormap');
% clim = [2,5];
% % significance level at t=3.7
% a = linspace(clim(1),clim(2),256);
% c(abs(a)<3.7,:) = repmat([1 1 1],nnz(abs(a)<3.7),1);
% 
% for ii = 2:length(fit_parameters)
%     
% [p,sind] = sort(pfit{ii});
% sigS = p'<(0.05./(length(channels):-1:1));
% % channels(sind(sigS)) % Show names of signficant channels
% %     grid on
% subplot(1,3,ii)
% T.trial{1} = Tfit{ii}; T.avg = Tfit{ii};% T.mask = sigS(sind)';
% cfg = [];
% cfg.channel = channels;
% cfg.layout = 'CTF275_helmet.mat';
% cfg.parameter = 'avg';
% cfg.interpolatenan = 'no';
% % cfg.maskparameter = 'avg';
% cfg.colormap   = c;
% cfg.zlim =clim;
% cfg.comment    = 'no';
% ft_topoplotER(cfg, T)
% 
% % titlename = fit_parameters{ii};
% % k = strfind(titlename,'_');
% % titlename(k) = ' ';
% title(titlename)%
% if ii == 2
%     title('E_t')
% elseif ii ==3
%     title('\eta_t')
% end
% end
% 
% subplot(1,3,1)
% caxis(clim)
% colorbar

%%

% Read data header from one subject to obtain sensot positions
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


% Declerad to use 2-tailed t-test with 2,000 randperms
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
M = cell2mat(M);
alpha = 0.05; 
alpha = alpha/2; % two-tailed
nparam = 2; % number of predictors for multiple comparisons
Mp = sort(M,'descend');
threshp = Mp(round(alpha/nparam*size(clusternull,1)));
Mn = sort(M,'ascend');
threshn = Mn(round(alpha/nparam*size(clusternull,1)));
%%
% figure(6); clf; set(gcf,'color','w','position',[176 748 1262 400])
for ii = 1:2
[tfced] = matlab_tfce_transform_MEG(Tfit{ii}',H,E,dh,neighbours);
[tfcen] = matlab_tfce_transform_MEG(-Tfit{ii}',H,E,dh,neighbours);
tfced = tfced - tfcen;
Tplot = Tfit{ii};
Tplot(tfced<threshp & tfced>threshn) = 0;

c = get(gca,'colormap');
if min(Tplot)==0
    clim = [2,4]; 
else
    clim = [-4,4]; 
end
a = linspace(clim(1),clim(2),256);
cthresh = min(Tplot(Tplot>0));
c(abs(a)<cthresh,:) = repmat([1 1 1],nnz(abs(a)<cthresh),1);
T = struct;
T.label = channels;
T.time{1} = 300;
T.sampleinfo = [1 1];
T.trial{1} =Tplot; T.avg = T.trial{1};% T.mask = sigS(sind)';

figure; clf; set(gcf,'color','w')

cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.parameter = 'avg';
cfg.interpolatenan = 'no';
cfg.colormap   = c;
cfg.zlim =clim;
cfg.comment    = 'no';
% subplot(1,2,ii)
ft_topoplotER(cfg, T)
colorbar

% saveas(gcf,sprintf('~/matlab/figures/powersens_gamma_1e4perms_confirm_%s.png',fit_parameters{ii}))
end

% saveas(gcf,sprintf('~/matlab/figures/powersens_gamma_1e4perms_confirm_twotailedcax.png'))