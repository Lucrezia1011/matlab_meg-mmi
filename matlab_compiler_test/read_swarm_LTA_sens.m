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

%%
load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
meg = zeros(npoints,nrois,size(Xv,1));
if ~exist('meg','var')
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
addpath('~/matlab/matlab_compiler_test/')
cd(data_path)
command_list = [];
for m = 1:length(fit_parameters)

    fit_parameter = fit_parameters{m};
    outpath = [data_path,'/lme_',fit_parameter,'/'];
    for ii = 1:length(param_list)
        filename = sprintf('ROI_%s.csv',param_list{ii});
        if ~exist(sprintf('%s%s',outpath,filename),'file')
            
            command_list{end+1} = {
                meg_data_name,latent_vars_name,param_list{ii},...
                num2str(npoints),fit_parameter,outpath};
        end
    end
    
end

if ~isempty(command_list)
parfor ii = 1:length(command_list)
    minputs = command_list{ii};
    mmi_LTA_trials_new(minputs{1},minputs{2},minputs{3},minputs{4},minputs{5},minputs{6});
    
end
end
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

clusteraal = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = Xfit{m};
%     if strncmp(freq,'evoked_',7)
        r = X.ROI;
%     else
%         [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
%     end
    n = 0;
    clusteraal{m} = zeros(length(param_list),npoints);

    for iir = 1:length(param_list)

        x = X(r==iir,:);
        LME = x.tStat;
        
        TFCE = tfce2d(LME);

        clusteraal{m}(iir,:) = TFCE';
        % Equivalent to 
%         [tfced] = matlab_tfce_transform(LME,2,0.5,4,0.1);   
%         [tfcedn] = matlab_tfce_transform(-LME,2,0.5,4,0.1);   
%         tfced = tfced - tfcedn;
    end
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement
close all
pv = 0.05;
plotopt = 't';

% Xv.RPE_abs = abs(Xv.RPE);

M = length(fit_parameters);
Xm = zeros(size(Xv,1),M);
for m = 1:M
    Xm(:,m) = Xv.(fit_parameters{m});
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
if isempty(clusternull)
    clusternull = cell2mat(clusteraal);
    M = nrois;
    C = corr(meg');
    lambda = eig(C);
    % Effective number of independent variables
    Meff = 1 + (M-1)*(1 - var(lambda)/M);
    alpha = 1 - (1 - pv)^(1/Meff);
  
end
alpha = alpha/length(fit_parameters);
snull = sort(abs(clusternull(:)));
Nnull = size(snull,1);
uplim =  snull(ceil(Nnull*(1-alpha)),1);
lowlim = -uplim;


% positive and negative separately
% snull = sort(clusternull);
% Nnull = size(snull,1);
% lowlim = snull(floor(Nnull*alpha),2);
% uplim =  snull(ceil(Nnull*(1-alpha)),1);

%%



% close all
dx = 1;
for m = [2:3,5:length(fit_parameters)]
    
    fit_parameter = fit_parameters{m};
    X = Xfit{m};
%     if strncmp(freq,'evoked_',7)
        r = X.ROI;
%     else
%         [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
%     end
    n = 0;

    for iir = 1:length(param_list)
     
        
%         clusternull = dlmread(['ROI_',param_list{iir},'_permute.txt']);
%         snull = sort(clusternull(:));
%         Nnull = size(snull,1);
%         lowlim = snull(floor(Nnull*pv/2));
%         uplim =  snull(ceil(Nnull*(1-pv/2)));     

        TFCE = clusteraal{m}(iir,:);
        
        x = X(r==iir,:);  
        if strcmp(plotopt,'t')
            megr = meg(iir,:);
            LME = x.tStat;           
        else
            megr = meg(iir,:);
            LME = x.Estimate; 
        end
        if any(TFCE>uplim)
            [pks,locs,w] = findpeaks(TFCE,1:npoints,'MinPeakHeight',uplim,'Annotate','extents');              
        else
            locs = [];
            w = [];
        end
        
        if any(TFCE<lowlim)
            [pksn,locsn,wn] = findpeaks(-TFCE,1:npoints,'MinPeakHeight',-lowlim,'Annotate','extents');
            locs = cat(2,locs,locsn);
            w = cat(2,w,wn);         
        end
        w = round(w);
        
        ind = cell(1,size(locs,2));
        for npks = 1:length(ind)
            ind{npks} = locs(npks)+ (-w(npks):w(npks));
        end
        
        if ~isempty(ind)
            
            if n == 0
                figure(m); clf
%                 set(gcf,'color','w','name',fit_parameter,'position',[346  464  1097  399])
                set(gcf,'color','w','name',fit_parameter,'position',[346  315  1055  548])
            end
            n = n+1;
            sb = subplot(3,5,n);
            yyaxis left
            plot(time,megr,'Linewidth',2); ylim(max(abs(ylim))*[-1 1])
            yyaxis right
            hold on
            plot(time,LME,'k')
            ylim(max(abs(ylim))*[-1 1])
            for npks = 1:length(ind)
                indpk = ind{npks};
                indpk(indpk<=0 | indpk > npoints) = [];
                plot(time(indpk),LME(indpk), 'r-','Linewidth',2)
            end

            if nrois == length(aal_labels)
                title(aal_labels{iir})
            elseif nrois == length(channels)
                title(channels{iir})
            end
%             xlabel('time (s)'); ylabel('tStat')
            grid on;
            if strcmp(freqb,'evoked')
                sb.XTick = -0.2:0.1:1;
                sb.XTickLabel([1,2,4,5,7,8,10,11,13]) = {''};          
                xlim([-.2,1]); 
            else
                sb.XTick = -3:1:3;
% %                 sb.XTickLabel([1,2,4,5,7,8,10,11,13]) = {''};
%                 sb.YTick = -6:2:6;
                xlim([-3,3]); 
            end
%             sb.XTickLabel = sb.XTick;
        end
            
    end
    if n ~= 0
%         saveas(gcf,sprintf('~/matlab/figures/%s_%s_final.tif',freq,fit_parameters{m}))
    end

end
