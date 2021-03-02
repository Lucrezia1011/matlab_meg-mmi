% read_swarm_LTA_final.m
% Lucrezia Liuzzi, last modified 2020/12/02
% Plot results of time point cluster analysis over ROIs 

clear all 
close all
clc
%% Read Swarm output

nrois = 116;
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/';


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

fit_parameters = Xv.Properties.VariableNames([3,5,7]);


addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
    
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

subs = unique(Xv.subject);
time = linspace(timew(1),timew(2),npoints);


%% Read source localized MEG data
if ~exist('meg','var')
    meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
    % reshape to: time x roi x trials
    meg = reshape(meg, npoints,nrois,size(meg,2));
    
    % baseline correct
    for r = 1:Xv.recording(end)
        ind = Xv.recording == r;
        base1 = mean(mean(meg(time<0,:,ind),3),1);
        base2 = meg(:,:,ind);
        meg(:,:,ind) = meg(:,:,ind) - base1;
    end
    
    % Plot high vs low surprise i.e. |RPE|
    S = std(mean(meg,3),0,1);
    [~,ind]= sort(S,'descend');
    %     c= corr(mean(meg(:,ind(1:30),:),3));
    %     aal_inds = {[1,4,10]; [2,3,7,8]; [5,6,9]}; % ROIs with largest variance over all trials
    %     aal_inds = {c(:,1)>0.9; c(:,2)>0.9; c(:,16)>0.9};
    
    figure; set(gcf,'color','w','position', [285   318   1031   808])
    for ii = 1:5
        subplot(2,5,ii); hold off
        %     inds = find(aal_inds{ii},1,'first');
        inds = ii;
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
        axis([-.2 1 -0.4 0.7])
        legend('|RPE|>4','|RPE|<2','diff','location','best')
        title(sprintf('Surprise effect, %s',aal_labels{ind(inds)}))
        grid on
    end
    
    % Plot positive vs negative RPE
    for ii = 1:5
        subplot(2,5,ii+5); hold off
        %     inds = find(aal_inds{ii},1,'first');
        inds = ii;
        ind1 = (Xv.RPE)>4;
        ind2 = (Xv.RPE)<-4;
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
        axis([-.2 1 -.4 .7])
        legend('RPE>4','RPE<-4','diff','location','best')
        title(sprintf('Value effect, %s',aal_labels{ind(inds)}))
        grid on
    end
%     saveas(gcf,'~/matlab/figures/LTA_AAL.png')
    
%     % Plot positive vs negative RPE, 1 ROI only
%     ii =40; % ROI
%     figure; clf
%     inds = ii;
%     ind1 = (Xv.RPE)>4;
%     ind2 = (Xv.RPE)<-4;
%     plot(time,mean(meg(:,ind(inds),ind1),3)); % high uncertainty
%     hold on
%     plot(time,mean(meg(:,ind(inds),ind2 ),3)); % low uncertainty
%     plot(time,mean(meg(:,ind(inds),ind1 ),3)-...
%         mean(meg(:,ind(inds),ind2 ),3),'k')
%     fill([time,fliplr(time)],...
%         [mean(meg(:,ind(inds),ind1 ),3)+std(meg(:,ind(inds),ind1),0,3)/sqrt(nnz(ind1));...
%         flipud(mean(meg(:,ind(inds),ind1 ),3)-std(meg(:,ind(inds),ind1 ),0,3)/sqrt(nnz(ind1)))],...
%         [0 0 1],'facealpha',.1,'edgecolor','none')
%     fill([time,fliplr(time)],...
%         [mean(meg(:,ind(inds),ind2 ),3)+std(meg(:,ind(inds),ind2),0,3)/sqrt(nnz(ind2));...
%         flipud(mean(meg(:,ind(inds),ind2 ),3)-std(meg(:,ind(inds),ind2 ),0,3)/sqrt(nnz(ind2)))],...
%         [1 0 0],'facealpha',.1,'edgecolor','none')
%     axis([-.2 1 -.4 .7])
%     legend('RPE>4','RPE<-4','diff','location','best')
%     title(sprintf('Value effect, %s',aal_labels{ind(inds)}))
%     grid on
    
    %%
    % Average over all trials for plotting final results
    meg = mean(meg,3);
    meg = meg';
    
end

% Plot average response from ROIs with largest variance
S = std(meg,0,2);
[~,ind]= sort(S,'descend');
figure; set(gcf,'color','w','position', [285   518   1031   408])

c = corr(meg(ind(1:16),:)');
subplot(131)
inds = [1,4,10];
plot(time,meg(ind(inds),:)');
legend(aal_labels{ind(inds)})
grid on; axis([-.2 1 -.3 0.6])
inds = [2,3,7,8];
subplot(132)
plot(time,meg(ind(inds),:)');
legend(aal_labels{ind(inds)})
grid on; axis([-.2 1 -.3 0.6])
inds = [5,6,9];
subplot(133)
plot(time,meg(ind(inds),:)');
legend(aal_labels{ind(inds)})
grid on; axis([-.2 1 -.3 0.6])


%% Check missing
% addpath('~/matlab/matlab_compiler_test/')
% cd(data_path)
% command_list = [];
% for m = 1:length(fit_parameters)
% 
%     fit_parameter = fit_parameters{m};
%     outpath = [data_path,freq,'/lme_',fit_parameter,'/'];
%     for ii = 1:length(param_list)
%         filename = sprintf('ROI_%s.csv',param_list{ii});
%         if ~exist(sprintf('%s%s/lme_%s/%s',data_path,freq,fit_parameter,filename),'file')
%             
%             command_list{end+1} = {
%                 meg_data_name,latent_vars_name,param_list{ii},...
%                 num2str(npoints),fit_parameter,outpath};
%         end
%     end
%     
% end

% if ~isempty(command_list)
% parfor ii = 1:length(command_list)
%     minputs = command_list{ii};
%     mmi_LTA_trials_new(minputs{1},minputs{2},minputs{3},minputs{4},minputs{5},minputs{6});
%     
% end
% end
%% Read statistics results from file
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
    clusteraal{m} = zeros(length(param_list),npoints);

    for iir = 1:length(param_list)

        x = X(r==iir,:);
        LME = x.tStat;
        
        TFCE = tfce2d(LME);
        % Equivalent to:
%         [tfced] = matlab_tfce_transform(LME,2,0.5,4,0.1);   
%         [tfcedn] = matlab_tfce_transform(-LME,2,0.5,4,0.1);   
%         tfced = tfced - tfcedn;
        clusteraal{m}(iir,:) = TFCE';
    end
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement
close all
pv = 0.05;
plotopt = 't';

% Only testing mood, RPE_LTA and RPE_sum
% M = 3;  %length(fit_parameters);
% Xm = zeros(size(Xv,1),M);
% for m = 1:M
%     Xm(:,m) = Xv.(fit_parameters{m});
% end
% C = corr(Xm);
% lambda = eig(C);
% % Effective number of independent variables
% Meff = 1 + (M-1)*(1 - var(lambda)/M);
% alpha = 1 - (1 - pv)^(1/Meff);


aal_confirm = [];

for m = 1:3
    
    aal_confirm(m).fit_parameter = fit_parameters{m};
    
    if m == 1
        aal_confirm(m).ROI = {'030'};
    else
        aal_confirm(m).ROI = {'068';'070'};
    end
end

% cannot combine over condition, because number of trials would be different  
% but can combine over predictor? 
ind = strfind(freq,'_');
freqb = freq(1:(ind-1));
clusternull = cell(length(fit_parameters),1);
condition = freq((ind+1):end);
m_ind = 0;
for m = 1:length(fit_parameters)
    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s_%s/lme_%s',data_path,freqb,condition,fit_parameter))
    
   
    % Each file records maximum and minimum TFCE values over all 116 ROIs,
    % correcting for multiple comparisons over all ROIs
    
    for d = 1:length(aal_confirm(m).ROI)
        dir_name = ['ROI',aal_confirm(m).ROI{d}(2:3),'_permute'];
        nullnames = dir(dir_name);
        nullnames(1:2)=[];
        clusternull2 = cell(length(nullnames),1);
        for n = 1:length(nullnames)
            if nullnames(n).bytes == 0
                delete([dir_name,'/',nullnames(n).name])
            else
                clusternull2{n} = dlmread([dir_name,'/',nullnames(n).name]);
            end
        end
        m_ind = m_ind +1;
        clusternull{m_ind} = cell2mat(clusternull2);
        %         clusternull{m} = cat(1,clusternull{m},cell2mat(clusternull2));

        alpha = 0.05/length(aal_confirm(m).ROI); % multiple comparisons
        snull = sort(clusternull{m_ind});
        Nnull = size(snull,1);
        aal_confirm(m).lowlim(d) = snull(floor(Nnull*alpha),2);
        aal_confirm(m).uplim(d) =  snull(ceil(Nnull*(1-alpha)),1);
    end
    
end

%% 
figure; subplot(311); 
histogram(clusternull{1})
title('Null distribution: mood')
subplot(312); 
histogram(clusternull{2})
hold on
histogram(clusternull{3})
title('Null distribution: RPE_LTA')

subplot(313); 
histogram(clusternull{4})
hold on
histogram(clusternull{5})
title('Null distribution: RPE_sum')
%%

% RPE and self-reported mood will be correlated with the variability in the
% evoked response to the gamble outcome (feedback), in right precuneus {68} and 
% paracentral lobule {70}(at ~500ms) for RPE and in the right insular cortex {30}
% (at ~400ms after feedback presentation) for mood. 

% close all
clc
dx = 1;
n=0;
for m = 1:3
    
    fit_parameter = fit_parameters{m};
    X = Xfit{m};
    r = X.ROI;
%     n = 0;
    
    for dd = 1:length(aal_confirm(m).ROI)
     
        iir =  str2double(aal_confirm(m).ROI{dd});
        TFCE = clusteraal{m}(iir,:);
        
        x = X(r==iir,:);  
        if strcmp(plotopt,'t')
            megr = meg(iir,:);
            LME = x.tStat;           
        else
            megr = meg(iir,:);
            LME = x.Estimate; 
        end
        
       
        %  Highlight peaks > TFCE null
            if any(TFCE>aal_confirm(m).uplim(dd))
                [pks,locs,w] = findpeaks(TFCE,1:npoints,'MinPeakHeight',aal_confirm(m).uplim(dd),'Annotate','extents');   
                [~,pksind] = max(pks);
                fprintf('Peaks of %s:\n%s at %.0fms,  T=%.2f,  p-value=%.1e\n',fit_parameter,...
                    aal_labels{iir}, time(locs(pksind))*1e3, LME(locs(pksind)),x.pValue(locs(pksind)) )
            else
                locs = [];
                w = [];
            end

            if any(TFCE<aal_confirm(m).lowlim(dd))
                [pksn,locsn,wn] = findpeaks(-TFCE,1:npoints,'MinPeakHeight',-aal_confirm(m).lowlim(dd),'Annotate','extents');
                fprintf('Peaks of %s:\n%s at %.0fms,  T=%.2f,  p-value=%.1e\n',fit_parameter,...
                        aal_labels{iir}, time(locsn)*1e3, LME(locsn),x.pValue(locsn) )
                locs = cat(2,locs,locsn);
                w = cat(2,w,wn);         
            end
            w = round(w);

            indc = cell(1,size(locs,2));
            for npks = 1:length(indc)
                indc{npks} = locs(npks)+ (-w(npks):w(npks));
            end

    %           Highlight timepoints > TFCE null
            indp = find(TFCE > aal_confirm(m).uplim(dd) | TFCE < aal_confirm(m).lowlim(dd));

            d = find(diff(indp)>1);
            d0 = 1;

            if nnz(indp) > 0
                ind = cell(1,nnz(d)+1);
                if isempty(d)
                    ind{1} =  indp;
                    fprintf('Peaks of %s:\n%s at %.0f-%.0fms\n\n',fit_parameter,...
                        aal_labels{iir}, time(indp(1))*1e3, time(indp(end))*1e3)
                else
                    for dii = 1:nnz(d)
                        ind{dii} = indp(d0:d(dii));
                        fprintf('Peaks of %s:\n%s at %.0f-%.0fms\n',fit_parameter,...
                            aal_labels{iir}, time(indp(d0))*1e3, time(indp(d(dii)))*1e3)

                        d0 = d(dii)+1;
                    end
                    ind{dii+1} = indp(d0:end);
                    fprintf('Peaks of %s:\n%s at %.0f-%.0fms\n\n',fit_parameter,...
                        aal_labels{iir}, time(indp(d0))*1e3, time(indp(end))*1e3)

                end
            else
                ind = [];
            end
        
        
        % Plot if there is a significant result or if we mention ROI in the
        % pre-registered hypotheses
            
            if n == 0
                figure(1); clf
                set(gcf,'color','w','name',fit_parameter,'position',[cd ~/mat])
%                 set(gcf,'color','w','name',fit_parameter,'position',[346  315  1055  190]) % 548
            end
            n = n+1;
            sb = subplot(1,5,n);
            yyaxis left
            plot(time,megr,'Linewidth',2); ylim(max(abs(ylim))*[-1 1])
            ylim([-0.55 0.55])
            sb.YTick = -.4:.2:.4;
            yyaxis right
            hold on
            % Plot whole cluster
            for npks = 1:length(indc)
                indpk = indc{npks};
                indpk(indpk<=0 | indpk > npoints) = [];
                plot(time(indpk),LME(indpk), '-','Linewidth',3 , 'color', [1 0.7 0.7])
            end
            
            plot(time,LME,'k-')
            ylim(max(abs(ylim))*[-1 1])
            
            % Plot peaks above threshold
            for npks = 1:length(ind)
                indpk = ind{npks};
                indpk(indpk<=0 | indpk > npoints) = [];
                plot(time(indpk),LME(indpk), 'r-','Linewidth',3)
            end
            ylim([-5.5 5.5])
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
                sb.YTick = -4:2:4;
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
saveas(gcf,sprintf('~/matlab/figures/%s_final_confirm.tif',freq))

