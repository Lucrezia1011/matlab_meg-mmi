clear all 
close all
clc
%% Read Swarm output

nrois = 116;
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/';

% nrois = 269;
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/';

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

freq = 'tfs_cue'; 

load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/24071_1')

if strcmp(freq(4:end),'cue')  
    time = tfscue.time;
    freqs = tfscue.freq;
elseif strcmp(freq(4:end),'choice')  
    time = tfschoice.time;
    freqs = tfschoice.freq;
else   
    time = tfsout.time;
    freqs = tfsout.freq;
end

clear Ycue Yout Ychoice ltvchoice ltvcue ltvout ...
    tfscue tfsout tfschoice ltvchoice_tfs ltvcue_tfs ltvout_tfs


meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
Xv = readtable([data_path,latent_vars_name],opts);

if contains(freq,'outcome')
    fit_parameters = Xv.Properties.VariableNames([3:9,11]);
else
    fit_parameters = Xv.Properties.VariableNames(3:9);
end

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
addpath ~/matlab
% if ~exist('atlas','var')    
%     atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
% end
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

subs = unique(Xv.subject);
meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
meg = reshape(meg, length(time),length(freqs),nrois,size(meg,2));
meg0 = meg;
for r = 1:18
    ind = Xv.recording == r;
%     base1 = mean(mean(meg(:,:,:,ind),4),1);
    base1 = mean(mean(meg(time<=0,:,:,ind),4),1);
    
%     meg(:,:,:,ind) = (meg(:,:,:,ind) - base1)./base1;
    meg(:,:,:,ind) = (meg(:,:,:,ind) - base1);
end


%%
ii = 69; 

cax = [-15 15];
figure(1);clf; set(gcf,'Name',aal_labels{ii})
subplot(221)
ind =  Xv.RPE<-1;
pcolor(time,freqs,mean(meg(:,:,ii,ind),4)')
shading flat; colorbar; caxis(cax)
title(sprintf('%.0f Negative Expectation & RPE<0',nnz(ind)))
subplot(222)
ind = Xv.E<2 & Xv.RPE > 2;
pcolor(time,freqs,mean(meg(:,:,ii,ind),4)')
shading flat; colorbar; caxis(cax)
title(sprintf('%.0f Negative Expectation & RPE>0',nnz(ind)))
subplot(223)
ind = Xv.E>0 & Xv.RPE >4;
pcolor(time,freqs,mean(meg(:,:,ii,ind ),4)')
shading flat; colorbar; caxis(cax)
title(sprintf('%.0f Positive Expectation & RPE > 2',nnz(ind)))
subplot(224)
ind = Xv.E>3 & Xv.RPE_abs<1;
pcolor(time,freqs,mean(meg(:,:,ii,ind),4)')
shading flat; colorbar; caxis(cax)
title(sprintf('%.0f Low risk gamble',nnz(ind)))

colormap jet

%% Check missing
addpath('~/matlab/matlab_compiler_test/')
cd(data_path)
command_list = [];
for m = 1:length(fit_parameters)

    fit_parameter = fit_parameters{m};
    
    for ii = 1:length(param_list)
        filename = sprintf('ROI_%s.csv',param_list{ii});
        if ~exist(sprintf('%s%s/lme_%s/%s',data_path,freq,fit_parameter,filename),'file')
            
            command_list{end+1} = {
                meg_data_name,latent_vars_name,param_list{ii},...
                num2str(npoints),fit_parameter,freq,data_path};
        end
    end
    
end

if ~isempty(command_list)
parfor ii = 1:length(command_list)
    minputs = command_list{ii};
    mmi_LTA_trials_new(minputs{1},minputs{2},minputs{3},minputs{4},minputs{5},minputs{6},minputs{7});
    
end
end
%% Read data from file
npoints =  length(time)*length(freqs);

Xfit = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = [];

    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s/lme_%s',data_path,freq,fit_parameter))
    
    for ii = 1:length(param_list)
        filename = sprintf('ROI_%s.csv',param_list{ii});
        opts = detectImportOptions(filename);
        x = readtable(filename,opts);
        x.index = x.index + (ii-1)*npoints; % adjust index
        X = cat(1,X,x);
    end
    Xfit{m} = X;
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement

clusteraal = cell(1,length(fit_parameters));

for m = 1:length(fit_parameters)
    X = Xfit{m};
    r = X.ROI;
    n = 0;
    clusteraal{m} = zeros(length(param_list),length(time), length(freqs));

    for iir = 1:length(param_list)
        
        x = X(r==iir,:);
        LME = x.tStat;
        LME = reshape(LME, length(time), length(freqs));       
        
        TFCE = tfce2d(LME);
        clusteraal{m}(iir,:,:) = TFCE;
    end
end

clusterallc = clusteraal;
Xfitc = Xfit;

%%
clusteraal = clusterallc;
Xfit = Xfitc;
for ii = 1:length(fit_parameters)
    clusteraal{ii}(91:end,:,:) = []; % Does not consider cerebellum
%     clusteraal{ii}(1:90,:,:) = [];
end
ind = [3,5,7];
clusteraal(ind)= []; fit_parameters(ind) =[];
Xfit(ind)= [];
%% Calculate Threshold free cluster enhancement

pv = 0.05;
plotopt = 't';

M = length(fit_parameters);
Xm = zeros(size(Xv,1),M);
for m = 1:M
    Xm(:,m) = eval(['Xv.',fit_parameters{m}]);
end
C = corr(Xm);
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

% cannot combine over condition, because number of trials would be different  
% but can combine over predictor? 
ind = strfind(freq,'_');
freqb = freq(1:(ind-1));
clusternull = cell(length(fit_parameters),1);
condition = freq((ind+1):end);
for m = 1:length(fit_parameters)
    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s_%s/lme_%s',data_path,freqb,condition,fit_parameter))
    if exist('ROI_permute2.txt','file')
        clusternull{m} = dlmread('ROI_permute2.txt');      
    end
    if exist('ROI_permute2','dir')
        nullnames = dir('ROI_permute2');
        nullnames(1:2)=[];
        clusternull2 = zeros(length(nullnames),2);
        for n = 1:length(nullnames)
            clusternull2(n,:) = dlmread(['ROI_permute2/',nullnames(n).name]);
        end
        clusternull{m} = cat(1,clusternull{m},clusternull2);
    end
    
end

clusternull = cell2mat(clusternull);

if isempty(clusternull)
    clusternull = cell2mat(clusteraal);
    alpha = alpha/58;
end
snull = sort(abs(clusternull(:)));
Nnull = size(snull,1);
uplim =  snull(ceil(Nnull*(1-alpha)),1);
lowlim = -uplim;

cc = colormap(jet);
cc(32:33,:) = [1 1 1; 1 1 1];
close all
for m = 1:length(fit_parameters)
    
    fit_parameter = fit_parameters{m};
    X = Xfit{m};
    r = X.ROI;
    n = 0;

    for iir = 1:90 %1:length(param_list)%1:90%

        TFCE = squeeze(clusteraal{m}(iir,:,:));
        
        x = X(r==iir,:);  
        LME = x.tStat;
        LME = reshape(LME, length(time), length(freqs));

        
        ind = abs(TFCE)>uplim;
        
        if nnz(ind)>0
            
            LMEplot = zeros(size(LME));
            zz= find(TFCE>uplim);
            S=regionprops(TFCE>uplim*.1,'PixelIdxList','PixelList');
            for z = 1:length(S)
                if ~isempty(intersect(zz,S(z).PixelIdxList))
                    LMEplot(S(z).PixelIdxList)  = LME(S(z).PixelIdxList);
                end
            end
            
            zz= find(TFCE<-uplim);
            S=regionprops(TFCE<-uplim*.2,'PixelIdxList','PixelList');
            for z = 1:length(S)
                if ~isempty(intersect(zz,S(z).PixelIdxList))
                    LMEplot(S(z).PixelIdxList)  = LME(S(z).PixelIdxList);
                end
            end
            
%             LME(abs(TFCE)<uplim) = 0;
            
            
            if n == 0
                figure(m); clf
%                 set(gcf,'color','w','name',fit_parameter,'position',[346  464  1097  399])
                set(gcf,'color','w','name',fit_parameter,'position',[346  315  1055  548])
                colormap(cc)
            end
            n = n+1;
            sb = subplot(3,5,n);
            pcolor(time,freqs,LMEplot'); shading flat;
            caxis([-3 3])
            if nrois == length(aal_labels)
                title(aal_labels{iir})
            elseif nrois == length(channels)
                title(channels{iir})
            end

        end
            
    end
    if n ~= 0
%         saveas(gcf,sprintf('~/matlab/figures/%s_%s.tif',freq,fit_parameters{m}))
    end

end

%%
iir = 88;
x = X(r==iir,:);
LME = x.tStat;
LME = reshape(LME, length(time), length(freqs));

% figure;
subplot(131)
pcolor(time,freqs,mean(meg(:,:,iir,:),4)'); shading flat;
caxis([-1 1])
subplot(132)
pcolor(time,freqs,mean(meg(:,:,iir,Xv.E_LTA>2.5 ),4)' - mean(meg(:,:,iir,Xv.E_LTA<2.5 ),4)'); shading flat;
caxis([-.5 .5])
subplot(133)
pcolor(time,freqs,LME'); shading flat;
caxis([-3 3])