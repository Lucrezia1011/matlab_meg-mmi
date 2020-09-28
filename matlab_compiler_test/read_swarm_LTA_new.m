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

freq = 'evoked_outcome'; 
% freq = 'theta_outcome'; 

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
fit_parameters = Xv.Properties.VariableNames([3:9]);


addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
    
% if ~exist('atlas','var')    
%     atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
% end
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

subs = unique(Xv.subject);
time = linspace(timew(1),timew(2),npoints);

%%

if ~exist('meg','var')
    meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
    meg = reshape(meg, npoints,nrois,size(meg,2));

% %%    
% 

    for r = 1:18
        ind = Xv.recording == r;
        base1 = mean(mean(meg(time<0,:,ind),3),1);   
        base2 = meg(:,:,ind);
        meg(:,:,ind) = (meg(:,:,ind) - base1)./std(base2(:));
    end
% 
% %%  
% % Cuneus 45 and paracentral lobule 69 have largest response to cue
%     
%     ii = 69; clf
% %     nneg = true(size(Xv.recording)); 
% %     npos = Xv.recording==6;
%     nneg = Xv.E_sum>14; 
%     npos = Xv.E_sum<14;
%     
%     subplot(211)
%     plot(time, mean(meg(:,ii,npos),3))
%     hold on
%     plot(time, mean(meg(:,ii,nneg),3)); title(aal_labels{ii})
%     fill([time,fliplr(time)],[mean(meg(:,ii,npos),3)+std(meg(:,ii,npos),0,3)/sqrt(nnz(npos));...
%         flipud(mean(meg(:,ii,npos),3)-std(meg(:,ii,npos),0,3)/sqrt(nnz(npos)))],...
%         [0 0 1],'facealpha',0.3,'edgecolor','none')
%     fill([time,fliplr(time)],[mean(meg(:,ii,nneg),3)+std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg));...
%         flipud(mean(meg(:,ii,nneg),3)-std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg)))],...
%         [1 0 0],'facealpha',0.3,'edgecolor','none');
%     %ylim([-0.7 0.7]); 
%     grid on; legend('posititve','negative')
%     pause(1)
%     
%     ii = ii+1;
%     subplot(212)
%      plot(time, mean(meg(:,ii,npos),3))
%     hold on
%     plot(time, mean(meg(:,ii,nneg),3)); title(aal_labels{ii})
%     fill([time,fliplr(time)],[mean(meg(:,ii,npos),3)+std(meg(:,ii,npos),0,3)/sqrt(nnz(npos));...
%         flipud(mean(meg(:,ii,npos),3)-std(meg(:,ii,npos),0,3)/sqrt(nnz(npos)))],...
%         [0 0 1],'facealpha',0.3,'edgecolor','none')
%     fill([time,fliplr(time)],[mean(meg(:,ii,nneg),3)+std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg));...
%         flipud(mean(meg(:,ii,nneg),3)-std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg)))],...
%         [1 0 0],'facealpha',0.3,'edgecolor','none');
%     %ylim([-0.7 0.7]); 
%     grid on; legend('positive','negative')
%     
% %%    

% megc = meg(:,[91:end],:);
% megc = meg(:,[43:56],:);
% figure; clf; set(gcf,'color','white')
% for ii = 1:13
%     subplot(3,5,ii)
%     megcc = megc(:,(ii-1)*2+(1:2),:);
%     plot(time, mean(megcc,3))
%     hold on
%     fill([time,fliplr(time)], [mean(megcc(:,1,:),3) + std(megcc(:,1,:),0,3)/sqrt(size(megcc,3)); ...
%         flipud(mean(megcc(:,1,:),3) - std(megcc(:,1,:),0,3)/sqrt(size(megcc,3)))],[0 0 1],'facealpha',0.2,'edgealpha',0)
%     fill([time,fliplr(time)], [mean(megcc(:,2,:),3) + std(megcc(:,2,:),0,3)/sqrt(size(megcc,3)); ...
%         flipud(mean(megcc(:,2,:),3) - std(megcc(:,2,:),0,3)/sqrt(size(megcc,3)))],[1 0 0],'facealpha',0.2,'edgealpha',0)
% %     title(aal_labels{(ii-1)*2+91})
%     title(aal_labels{(ii-1)*2+43})
%     xlim([-0.5, 0.8]); grid on; ylim([-0.3 0.3])
% end

    meg = mean(meg,3);
    meg = meg';

end

% figure; set(gcf,'color','white')
% subplot(121)
% plot(time,meg(1,:),'LineWidth',2)
% hold on
% plot(time,meg([19,77,47],:))
% plot(time,mean(meg(98:2:100,:),1),'LineWidth',2)
% xlim([-.5, 0.8]); ylim([-.42 .42]); grid on
% title('Left cerebrum and right cerebellum')
% xlabel('Time (s)'); ylabel('Evoked response to choice (z-score)')
% plot([-0.022 -0.022],[-.5 .5],'-','color',[0.5 0.5 0.5])
% plot([0.088 0.088],[-.5 .5],'-','color',[0.5 0.5 0.5])
% plot([0.168 0.168],[-.5 .5],'-','color',[0.5 0.5 0.5])
% hold off
% 
% subplot(122)
% plot(time,meg(2,:),'LineWidth',2)
% hold on
% plot(time,meg([20,78,48],:))
% plot(time,mean(meg(97:2:99,:),1),'LineWidth',2)
% xlim([-.5, 0.8]); ylim([-.42 .42]); grid on
% legend('Precentral C.','Supp.Motor.C.','Thalamus','Lingual Gy.','Cerebellum','-22ms','88ms','168ms')
% title('Right cerebrum and left cerebellum')
% xlabel('Time (s)')
% plot([-0.022 -0.022],[-.5 .5],'-','color',[0.5 0.5 0.5])
% plot([0.088 0.088],[-.5 .5],'-','color',[0.5 0.5 0.5])
% plot([0.168 0.168],[-.5 .5],'-','color',[0.5 0.5 0.5])

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
             
    end
end

%%
fit_parameters([3,5,7]) = [];
clusteraal([3,5,7]) = [];
Xfit([3,5,7])=[];
%% Calculate Threshold free cluster enhancement
close all
pv = 0.05;
plotopt = 't';

Xv.RPE_abs = abs(Xv.RPE);

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

% clusternull = cell(1,length(param_list));
% 
% for iir = 1:length(param_list)   
%     clusternull{iir} = dlmread(['ROI_',param_list{iir},'_permute.txt']);
% end
% clusternull  = cell2mat(clusternull');

% Consider 2-tailed distribution
% snull = sort(clusternull);
% Nnull = size(snull,1);
% lowlim = snull(floor(Nnull*alpha),2);
% uplim =  snull(ceil(Nnull*(1-alpha)),1);

% Consider 1-tailed distribution
% snull = sort(abs(clusternull(:)));
% Nnull = size(snull,1);
% uplim =  snull(ceil(Nnull*(1-alpha)),1);
% lowlim = -uplim;

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
    alpha = 0.05/116;
   
end
snull = sort(abs(clusternull(:)));
Nnull = size(snull,1);
uplim =  snull(ceil(Nnull*(1-alpha)),1);
lowlim = -uplim;
% 
%%
% snull = sort(clusternull);
% Nnull = size(snull,1);
% lowlim = snull(floor(Nnull*alpha),2);
% uplim =  snull(ceil(Nnull*(1-alpha)),1);
close all
dx = 1;
for m = [1:3,5]%1:length(fit_parameters)
    
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
        saveas(gcf,sprintf('~/matlab/figures/%s_%s_new.tif',freq,fit_parameters{m}))
    end

end
