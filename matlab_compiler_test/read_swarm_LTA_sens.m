clear all 
close all
clc
%% Read Swarm output

nrois = 266;


param_list = cell(nrois,1);
for nn = 1:nrois
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

freq = 'evoked_outcome'; 
data_path = ['/data/MBDU/MEG_MMI3/results/mmiTrial_sens/',freq,'/'];


if strcmp(freq,'evoked_choice')
    npoints = 360;
    timew = [-.5,.7];
else
    npoints = 360;
    timew = [-.2,1];
end
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
Xv = readtable([data_path,latent_vars_name],opts);

% fit_parameters = Xv.Properties.VariableNames([3,5,7:9,11]);
fit_parameters = Xv.Properties.VariableNames([3,5,7:10]);


addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
    
% if ~exist('atlas','var')    
%     atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
% end
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

subs = unique(Xv.subject);
time = linspace(timew(1),timew(2),npoints);

load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
%%
if ~exist('meg','var')
    meg = zeros(npoints,nrois,size(Xv,1));
    for nn = 1:nrois
        meg(:,nn,:) = dlmread(sprintf('%smeg_trials/sens_%s.txt',data_path,param_list{nn}));
    end
    
    S = std(mean(meg,3),0,1);
    [~,ind]= sort(S,'descend');
    c= corr(mean(meg(:,ind(1:30),:),3));
    [coeff,score]=pca(mean(meg,3));
    t = .92;
    aal_inds = {c(:,1)>t; c(:,2)>t; c(:,9)>t; c(:,5)>.98; c(:,29)>.98}; % ROIs with largest variance over all trials
    figure; set(gcf,'color','w','position', [285   318   1031   808])
    for ii = 1:5
    subplot(2,5,ii); hold off
    inds = find(aal_inds{ii},1,'first');
    ind1 = abs(Xv.RPE)>4;
    ind2 = abs(Xv.RPE)<2;
    plot(time,mean(meg(:,ind(inds),ind1),3)); % high uncertainty
    hold on
    plot(time,mean(meg(:,ind(inds),ind2 ),3)); % low uncertainty
    plot(time,mean(meg(:,ind(inds),ind1 ),3)-...
        mean(meg(:,ind(inds),ind2 ),3),'k')
    fill([time,fliplr(time)],...
        [mean(meg(:,ind(inds),ind1 ),3)+std(meg(:,ind(inds),ind1),0,3)/sqrt(nnz(ind1));...
        flipud(mean(meg(:,ind(inds),ind1 ),3)-std(meg(:,ind(inds),ind1 ),0,3)/sqrt(nnz(ind1)))],...
        [0 0 1],'facealpha',.1,'edgecolor','none')
    fill([time,fliplr(time)],...
        [mean(meg(:,ind(inds),ind2 ),3)+std(meg(:,ind(inds),ind2),0,3)/sqrt(nnz(ind2));...
        flipud(mean(meg(:,ind(inds),ind2 ),3)-std(meg(:,ind(inds),ind2 ),0,3)/sqrt(nnz(ind2)))],...
        [1 0 0],'facealpha',.1,'edgecolor','none')
    axis([-.2 1 -12e-14 12e-14])
    legend('|RPE|>4','|RPE|<2','diff','location','best')
    title(sprintf('Surprise effect, %s',channels{ind(inds)}))
    grid on
    end
   
  
    for ii = 1:5
    subplot(2,5,ii+5); hold off
    inds = find(aal_inds{ii},1,'first');
    ind1 = (Xv.RPE)>4;
    ind2 = (Xv.RPE)<-1;
    plot(time,mean(meg(:,ind(inds),ind1),3)); % high uncertainty
    hold on
    plot(time,mean(meg(:,ind(inds),ind2 ),3)); % low uncertainty
    plot(time,mean(meg(:,ind(inds),ind1 ),3)-...
        mean(meg(:,ind(inds),ind2 ),3),'k')
    fill([time,fliplr(time)],...
        [mean(meg(:,ind(inds),ind1 ),3)+std(meg(:,ind(inds),ind1),0,3)/sqrt(nnz(ind1));...
        flipud(mean(meg(:,ind(inds),ind1 ),3)-std(meg(:,ind(inds),ind1 ),0,3)/sqrt(nnz(ind1)))],...
        [0 0 1],'facealpha',.1,'edgecolor','none')
    fill([time,fliplr(time)],...
        [mean(meg(:,ind(inds),ind2 ),3)+std(meg(:,ind(inds),ind2),0,3)/sqrt(nnz(ind2));...
        flipud(mean(meg(:,ind(inds),ind2 ),3)-std(meg(:,ind(inds),ind2 ),0,3)/sqrt(nnz(ind2)))],...
        [1 0 0],'facealpha',.1,'edgecolor','none')
    axis([-.2 1 -12e-14 12e-14])
    legend('RPE>4','RPE<-1','diff','location','best')
    title(sprintf('Value effect, %s',channels{ind(inds)}))
    grid on
    end
%     saveas(gcf,'~/matlab/figures/LTAsens.png')
    
% %%    
% 

%     for r = 1:Xv.recording(end)
%         ind = Xv.recording == r;
%         base1 = mean(mean(meg(time<0,:,ind),3),1);   
%         base2 = meg(:,:,ind);
%         meg(:,:,ind) = meg(:,:,ind) - base1;
%     end

    meg = mean(meg,3);
    meg = meg';

end

cfg = [];
cfg.channel      = channels;
cfg.ylim         = [-1 1]*0.5e-13;
cfg.xlim         = [-0.2 1];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF275_helmet.mat';
cfg.interactive  = 'no';
meg_plot = [];
meg_plot.label = channels;
meg_plot.time{1} = time;
meg_plot.trial{1} = meg;
wind_size =  [491  66  1030   821];
figure(1); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotER(cfg, meg_plot); 

%% Check missing
% addpath('~/matlab/matlab_compiler_test/')
% cd(data_path)
% command_list = [];
% for m = 1:length(fit_parameters)
% 
%     fit_parameter = fit_parameters{m};
%     outpath = [data_path,'/lme_',fit_parameter,'/'];
%     for ii = 1:length(param_list)
%         filename = sprintf('ROI_%s.csv',param_list{ii});
%         if ~exist(sprintf('%s%s',outpath,filename),'file')
%             
%             command_list{end+1} = {
%                 meg_data_name,latent_vars_name,param_list{ii},...
%                 num2str(npoints),fit_parameter,outpath};
%         end
%     end
%     
% end
% 
% if ~isempty(command_list)
% parfor ii = 1:length(command_list)
%     minputs = command_list{ii};
%     mmi_LTA_trials_new(minputs{1},minputs{2},minputs{3},minputs{4},minputs{5},minputs{6});
%     
% end
% end
%% Read data from file
Xfit = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = [];

    fit_parameter = fit_parameters{m};
    cd(sprintf('%s/lme_%s',data_path,fit_parameter))
    
    for ii = 1:length(param_list)
        filename = sprintf('ROI_%s.csv',param_list{ii});
        opts = detectImportOptions(filename);
        x = readtable(filename,opts);
        x.index = x.index + (ii-1)*npoints; % adjust index
        x.ROI = repmat(ii,npoints,1);
        X = cat(1,X,x);
    end
    Xfit{m} = X;
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement
% close all
pv = 0.05;
% 
% M = 3;%length(fit_parameters);
% Xm = zeros(size(Xv,1),M);
% for m = 1:M
%     Xm(:,m) = Xv.(fit_parameters{m});
% end
% C = corr(Xm);
% lambda = eig(C);
% % Effective number of independent variables
% Meff = 1 + (M-1)*(1 - var(lambda)/M);
% alpha = 1 - (1 - pv)^(1/Meff);
% cannot combine over condition, because number of trials would be different  
% but can combine over predictor? 
ind = strfind(freq,'_');
freqb = freq(1:(ind-1));
clusternull = cell(length(fit_parameters),1);
condition = freq((ind+1):end);
for m = 1:length(fit_parameters)
    fit_parameter = fit_parameters{m};
    cd(sprintf('%s/lme_%s',data_path,fit_parameter))
    if exist('ROI_permute.txt','file')
        clusternull{m} = dlmread('ROI_permute.txt');      
    end
    if exist('ROI_permute','dir')
        nullnames = dir('ROI_permute');
        nullnames(1:2)=[];
        clusternull2 = cell(length(nullnames),1);
        for n = 1:length(nullnames)
            if nullnames(n).bytes == 0
                delete(['ROI_permute/',nullnames(n).name])
            else
                clusternull2{n} = dlmread(['ROI_permute/',nullnames(n).name]);
            end
        end
        clusternull{m} = cat(1,clusternull{m},cell2mat(clusternull2));
    end
    
end

clusternull = cell2mat(clusternull);
alpha = 0.01; 
alpha = alpha/2; % two-tailed
nparam = 3; % number of predictors for multiple comparisons
M = clusternull(:);
Mp = sort(M,'descend');
threshp = Mp(round(alpha/nparam*size(clusternull,1)));
Mn = sort(M,'ascend');
threshn = Mn(round(alpha/nparam*size(clusternull,1)));

%% Calculate Threshold free cluster enhancement
% Read data header from one subject to obtain sensot positions
cd(['/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/'])
hdr = ft_read_header('sub-24071_task-mmi3_run-1_meg.ds');

cfg_neighb        = [];
cfg_neighb.method = 'distance';%'distance';
cfg_neighb.channel = channels;
cfg_neighb.template = 'CTF275_neighb.mat';
cfg_neighb.neighbourdist    = 3.5;
neighbours        = ft_prepare_neighbours(cfg_neighb, hdr);
for n = 1:length(neighbours)
    [~,~,iB] = intersect(neighbours(n).neighblabel, channels );
    neighbours(n).neighbnum =iB;
end

chPos = hdr.grad.chanpos(strcmp(hdr.grad.chantype,'meggrad'),:);
figure; plot3(chPos(:,1),chPos(:,2),chPos(:,3),'.')
for ii = 1:length(chPos)
    ind = false(length(chPos),1);
    ind(ii) = true;
    m(ii) = min(sqrt(sum((chPos(ind,:) - chPos(~ind,:)).^2,2)));
end

hold on
plot3(chPos(:,1),chPos(:,2),chPos(:,3),'o')

chPosr = round((chPos - min(chPos))/1.8)+1;
img = zeros([14,12,11]);
for ii = 1:length(chPos)
    img(chPosr(ii,1),chPosr(ii,2),chPosr(ii,3)) = 1;
end

% Declerad to use 2-tailed t-test with 2,000 randperms
E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1;

%%
TFCE = cell(1,length(fit_parameters));
for m =1:3%length(fit_parameters)
    X = Xfit{m};
    T = reshape(X.tStat,[360,length(channels)])';
    % Add edges inter channels between time??
    [tfced] = matlab_tfce_transform_MEGtime(T,H,E,dh,neighbours);
    [tfcen] = matlab_tfce_transform_MEGtime(-T,H,E,dh,neighbours);
    tfced = tfced - tfcen;
    TFCE{m} = tfced;
end

%% Plots
for m =1:3
fit_parameter = fit_parameters{m};
X = Xfit{m};
T = reshape(X.tStat,[360,length(channels)])';

tfced = TFCE{m};
% 
% 
% figure; subplot(121);imagesc(T)
% colorbar; caxis([-4.5 4.5])
% subplot(122)
% imagesc(tfced); colorbar; caxis([-400 400])



Tplot = T;
Tplot(tfced<threshp & tfced>threshn) = 0;


c = get(gca,'colormap');
clim = [-5,5];
a = linspace(clim(1),clim(2),256);
cthresh = min(Tplot(Tplot>0));
c(abs(a)<cthresh,:) = repmat([1 1 1],nnz(abs(a)<cthresh),1);
T = struct;
T.label = channels;
T.time{1} = linspace(-.2,1,360);
T.sampleinfo = [1 360];
T.trial{1} =Tplot; T.avg = Tplot;% T.mask = sigS(sind)';

figure(m+1); 
clf
set(gcf,'position', [163   60  1290  922],'color','w','name',fit_parameter)
ii = 0;

tstep = 0.02;
for tt = .15:tstep:0.7
    ii = ii +1;
    cfg = [];
    cfg.channel = channels;
    cfg.layout = 'CTF275_helmet.mat';
    cfg.parameter = 'avg';
    cfg.interpolatenan = 'no';
    cfg.colormap   = c;
    cfg.zlim =clim;
    cfg.xlim = tt+[0 tstep];
    cfg.comment    = 'no';
    subplot(4,7,ii)
    ft_topoplotER(cfg, T)
%     colorbar
    title(sprintf('%.0f-%.0fms',tt*1e3,(tt+tstep)*1e3))
end

% saveas(gcf,sprintf('~/matlab/figures/%s_%s_sens.tif',freq,fit_parameter))

end

%%

mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
% sourcemodel.coordsys = 'mni';
load([ftpath,'headmodel/standard_mri.mat'])

sub = '24071';% '23490';%'24071';
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'];
cd(data_path)
cfg = [];
cfg.channel  = channels;
cfg.dataset = ['sub-',sub,'_task-mmi3_run-1_meg.ds'];
cfg.trials  = {'1'};
data = ft_preprocessing(cfg);

mri_name = [data_path(1:end-4),'anat/sub-',sub,'_acq-mprage_T1w.nii'];
if ~exist(mri_name,'file')
    mri_name = [mri_name,'.gz'];
end
fids_name =  ['sub-',sub,'_fiducials.tag'];
mri = fids2ctf(mri_name,fids_name,0);

cd(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/sub-',sub,'_task-mmi3_run-1_meg'])
load('headmodel.mat'); %vol
h = load('leadfields_5mm.mat');

m = 1;
% X = Xfit{m};
% Tplot = reshape(X.tStat,[360,length(channels)])';
Tplot= TFCE{m};
T = struct;
T.label = channels;
T.time{1} = linspace(-.2,1,360);
T.sampleinfo = [1 1];
timesel = [0.6500 0.700];
T.time{1} = mean(timesel);
time = linspace(-.2,1,360);
[~,ii]= min(abs(time-timesel'),[],2);
tt = ind2sub([360,2],ii);
Tplot = mean(Tplot(:,tt(1):tt(2)),2);
Tplot = Tplot*1e-15;
T.trial{1} =Tplot; T.avg = Tplot;
T.grad = data.grad;
T.hdr  =data.hdr;


cfg = [];
cfg.numdipoles  = 1; % number, default is 1
cfg.symmetry    = [];%'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
cfg.channel     = channels;
cfg.gridsearch  = 'yes';
cfg.resolution  = 8; % in mm
% cfg.sourcemodel = h.grid;
cfg.headmodel  = vol;
[source] = ft_dipolefitting(cfg, T);
pos = source.dip.pos;
% 
% 
% figure
% hold on
% 
% ft_plot_dipole(source.dip.pos(1,:), mean(source.dip.mom(1:3,:),2), 'color', 'r')
% pos = source.dip.pos;
% mri_resliced_cm = mri;
% ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1);
% ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
% ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);
% view(60,30)
% ft_plot_crosshair(pos, 'color', [1 1 1]/2);
% 
% axis tight
% axis off
% 
% % print -dpng natmeg_dip_sourcedif.png

cfg = [];
cfg.location = pos(1,:);
ft_sourceplot(cfg, mri)