
clear all 
close all
clc
%% Read Swarm output


nrois = 116;
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_old/';

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

freq = 'delta_outcome'; 
% freq = 'theta_outcome'; 

if strncmp(freq,'evoked_',7)
    npoints = 360;
    timew = [-.2,1];
else
    npoints = 360;
    timew = [-3,3];
end
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


% if strcmp(freq,'evoked')
%     npoints = 300;
%     meg_data_name = 'meg_trials.txt';
%     timew = [-.5,1];
%     latent_vars_name = 'latent_vars.csv';
% elseif strncmp(freq,'evoked_',7)
%     npoints = 360;
%     timew = [-.2,1];
%     meg_data_name = ['meg_trials_',freq,'.txt'];
%     latent_vars_name = ['latent_vars_',freq,'.csv'];
% 
% elseif strcmp(freq,'alpha')
%     npoints = 600;
%     meg_data_name = sprintf('meg_trials_%s.txt',freq);
%     timew = [-3,3];
%     latent_vars_name = 'latent_vars_3s.csv';
% 
% else
%     npoints = 150;
%     meg_data_name = sprintf('meg_trials_%s.txt',freq);
%     timew = [-.5,1];
%     latent_vars_name = 'latent_vars.csv';
% 
% end
%    

opts = detectImportOptions([data_path,latent_vars_name]);
Xv = readtable([data_path,latent_vars_name],opts);

fit_parameters = Xv.Properties.VariableNames(3:end);
fit_parameters([6,7])= [];


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
if ~exist('channels','var')
    % Find common channels
    sn = 2;
    sub = subn(sn,:);
    data_paths = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_paths)
    data_name = [sub,'MMI_mmi3_proc.ds'];
    h = ft_read_sens(data_name,'senstype','MEG');
    label2 = h.label(strncmp(h.label,'M',1));

    sn = 4;
    sub = subn(sn,:);
    data_paths = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_paths)
    data_name = [sub,'MMI_mmi3_proc.ds'];
    h = ft_read_sens(data_name,'senstype','MEG');
    label4 = h.label(strncmp(h.label,'M',1));

    channels = intersect(label2,label4);
end

subs = unique(Xv.subject);
if ~exist('meg','var')
    meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
    % find correlation at the subject level
    meg = reshape(meg, npoints,nrois,size(meg,2));
    Cs = zeros(nrois,nrois,length(subs));
    for s = 1:length(subs)
        megs = meg(:,:,Xv.subject==subs(s));
        Cs(:,:,s) = corr(mean(megs,3));
    end
    
    meg = mean(meg,3);
    meg = meg';
    
%     figure; set(gcf,'color','w')
%     subplot(221); 
%     C = mean((Cs),3);
%     lambda = eig(C); Meff = 1 + (nrois-1)*(1 - var(lambda)/nrois);
%     imagesc(C);  colorbar; caxis([-1,1]);
%     xlabel('ROI'); ylabel('ROI'); 
%     title(sprintf('Corr of average evoked responses per subject\nVar_l = %.2f , M_{eff} = %.1f',var(lambda),Meff));
% 
%     subplot(222); 
%     C = mean(abs(Cs),3);
%     lambda = eig(C); Meff = 1 + (nrois-1)*(1 - var(lambda)/nrois);
%     imagesc(C);  colorbar; caxis([0,1]);
%     xlabel('ROI'); ylabel('ROI'); 
%     title(sprintf('Abs Corr of average evoked responses per subject\nVar_l = %.2f , M_{eff} = %.1f',var(lambda),Meff));
%     
%     subplot(223); 
%     C = corr(meg');
%     lambda = eig(C); Meff = 1 + (nrois-1)*(1 - var(lambda)/nrois);
%     imagesc(C);  colorbar; caxis([-1,1]);
%     xlabel('ROI'); ylabel('ROI'); 
%     title(sprintf('Corr of average evoked responses \nVar_l = %.2f , M_{eff} = %.1f',var(lambda),Meff));
% 
%     subplot(224); 
%     C = abs(corr(meg'));
%     lambda = eig(C); Meff = 1 + (nrois-1)*(1 - var(lambda)/nrois);
%     imagesc(C); colorbar; caxis([0,1]);
%     xlabel('ROI'); ylabel('ROI'); 
%     title(sprintf('Abs Corr of average evoked responses \nVar_l = %.2f , M_{eff} = %.1f',var(lambda),Meff));
    
%     meg = mean(meg,2);
%     if strncmp(freq,'evoked_',7) || regexp(freq,'outcome')
%         meg = reshape(meg, npoints,nrois);
%         meg = meg';
%     else
%         meg = reshape(meg, nrois,npoints);
%     end
end
time = linspace(timew(1),timew(2),npoints);

%% Check for missing points
% 
% n = 0;
% runcompiled = 'run_mmi_LTA_trials.sh';               
% 
% command_list = [];   
% for m = 1:length(fit_parameters)
%     fit_parameter = fit_parameters{m};
%     cd(sprintf('%s/%s/lmixmodel_%s',data_path,freq,fit_parameter))
%     for ii = 1:length(param_list)
%         filename = sprintf('ROI_%s.csv',param_list{ii});
%         
%         if ~exist(filename,'file')
%             n = n +1;
%             command_list = [command_list ...
%                 sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID; '...
%                 '  test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz '...
%                 '  && ~/matlab/matlab_compiler_test/%s '...
%                 ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s %s\n'],runcompiled,...
%                 meg_data_name,latent_vars_name,param_list{ii},num2str(npoints),fit_parameter,freq,data_path)];
%             
%         end
%         
%     end
% end
% 
% % identify missing parameters: run as normal function
% cd ~/matlab/matlab_compiler_test/
% fit_parameter = 'RPE';
% param_list = {'022';'049';'053';'058';'070';'072';'073'};
% 
% parfor ii = 1:n
%     mmi_LTA_trials(meg_data_name,latent_vars_name,param_list{ii},num2str(npoints),fit_parameter,freq,data_path)
% end

% run as swarm
% if ~isempty(command_list)
%     cd ~/matlab/matlab_compiler_test
% 
%     file_handle = fopen(sprintf('mmi_LTA_trials_missed.swarm'),'w+');
%     fprintf(file_handle,command_list);
%     fclose(file_handle);
%     
%     emailnote = '"--mail-type=FAIL,END"';
%     % need to include lscratch! see matlab biowulf page
%     mem = '2';  % gigabytes
%     threads = '4'; % number of threads
%     bundles = '10'; % limits number of jobs running at the same time
%     logfolder = '~/matlab/matlab_compiler_test/swarm_logs';
%     jobid = evalc(sprintf('!swarm --job-name lmixmiss --gres lscratch:20 -g %s -t %s -b %s --time 00:30:00 --logdir %s -f mmi_LTA_trials_missed.swarm --sbatch %s',...
%         mem,threads,bundles,logfolder,emailnote));
%     
% end

%% Read data from file
Xfit = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = [];

    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s/lmixmodel_%s',data_path,freq,fit_parameter))
    
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
E = 0.5;
H = 2;

clusteraal = cell(1,length(fit_parameters));
peakaream = cell(1,length(fit_parameters));
dx = 1;
for m = 1:length(fit_parameters)
    X = Xfit{m};
%     if strncmp(freq,'evoked_',7)
        r = X.ROI;
%     else
%         [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
%     end
    n = 0;
    clusteraal{m} = zeros(length(param_list),npoints);
    peakarea = cell(1,length(param_list));

    for iir = 1:length(param_list)
        
        TFCE = zeros(1,npoints);

        x = X(r==iir,:);
        LME = x.tStat;
        
        indp = LME>0;
        ii = find(diff(indp));
        if indp(1) == 0 && indp(end) == 0
            a = cell(1,length(ii)/2);
            for jj = 1:length(ii)/2
                a{jj} = ii(jj*2-1)+1:ii(jj*2);
            end
        elseif indp(1) == 1 && indp(end) == 0
            a = cell(1,ceil(length(ii)/2));
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
        elseif indp(1) == 0 && indp(end) == 1
            a = cell(1,ceil(length(ii)/2));
            for jj = 1:floor(length(ii)/2)
                a{jj} = ii(jj*2-1)+1:ii(jj*2);
            end
            a{end} = ii(end)+1:npoints;
        else
            a = cell(1,length(ii)/2+1);
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)-1
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
            a{end} = ii(end)+1:npoints;
        end
        
        Ap = zeros(1,length(a));
        % convert into area
        for ii = 1:length(a)
            h = LME(a{ii});
            
            dh = 0.1;
            for p = 1:length(h)
                A = 0;
                hp = dh;
                while hp < h(p)
                    e = hp <= h;
                    de = diff(e);
                    de(end+1,:) = 0;
                    e = e & ~de;
                    for t = 1:length(e)
                        if e(t) == 0 && t<p
                            e(1:t) = 0;
                        elseif e(t) == 0 && t>p
                            e(t:end) = 0;
                        end
                    end
                    eE = (sum(e)*dx)^E;
                    
                    A = A + eE*(hp^H)*dh;
                    hp = hp + dh;
                end
                TFCE(p+a{ii}(1)-1) = A;
            end
            Ap(ii) = sum(h);
        end
        
        
        indp = LME<0;
        ii = find(diff(indp));
        if indp(1) == 0 && indp(end) == 0
            a = cell(1,length(ii)/2);
            for jj = 1:length(ii)/2
                a{jj} = ii(jj*2-1)+1:ii(jj*2);
            end
        elseif indp(1) == 1 && indp(end) == 0
            a = cell(1,ceil(length(ii)/2));
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
        elseif indp(1) == 0 && indp(end) == 1
            a = cell(1,ceil(length(ii)/2));
            for jj = 1:floor(length(ii)/2)
                a{jj} = ii(jj*2-1)+1:ii(jj*2);
            end
            a{end} = ii(end)+1:npoints;
        else
            a = cell(1,length(ii)/2+1);
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)-1
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
            a{end} = ii(end)+1:npoints;
        end
        
        An = zeros(1,length(a));
        % convert into area
        for ii = 1:length(a)
            h = -LME(a{ii});
            
            dh = 0.1;
            for p = 1:length(h)
                A = 0;
                hp = dh;
                while hp <= h(p)
                    e =  hp <= (h+dh/2);
                    de = diff(e);
                    de(end+1,:) = 0;
                    e = e & ~de;
                    for t = 1:length(e)
                        if e(t) == 0 && t<p
                            e(1:t) = 0;
                        elseif e(t) == 0 && t>p
                            e(t+1:end) = 0;
                        end
                    end
                    
                    %                 hold on
                    %                 plot(find(e),ones(1,nnz(e))*hp)
                    eE = (sum(e)*dx)^E;
                    
                    A = A + eE*(hp^H)*dh;
                    hp = hp + dh;
                end
                TFCE(p+a{ii}(1)-1) = -A;
            end
            An(ii) = sum(h);
        end
        
        clusteraal{m}(iir,:) = TFCE;
        peakarea{iir} = cat(2,An,Ap);
        
    end
   peakaream{m} = cell2mat(peakarea);
end

%% Plots with TFCE
% close all
% 
% pv = 0.01/5;
% sa = sort(cell2mat(peakaream));
% Athresh = sa(ceil((1-pv)*size(sa,2)));
% 
% stfce = cell2mat(clusteraal);
% stfce = sort(abs(stfce(:)));
% l = size(stfce,1);
% tresh = stfce(ceil((1-pv/2)*l));
% 
% for m = 1:length(fit_parameters)
%     fit_parameter = fit_parameters{m};
%     X = Xfit{m};
% %     if strncmp(freq,'evoked_',7) 
%         r = X.ROI;
% %     else
% %         [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
% %     end
%     n = 0;
%     
%     cd(sprintf('%s%s/lmixmodel_%s',data_path,freq,fit_parameter))
%     for iir = 1:length(param_list)
%         
%         x = X(r==iir,:);
%         LME = x.tStat;
% 
%         TFCE = clusteraal{m}(iir,:);
%         ind = find(abs(TFCE)>tresh );
%         if any(ind) %length(ind)>10  
%             
%             if n == 0
%                 figure(m); clf
%                 set(gcf,'color','w','name',fit_parameter,'position',[637   204   859   661])
% %                 set(gcf,'color','w','name',fit_parameter,'position',[1000 300 1350 950])
%             end
%             n = n+1;
%             sb = subplot(6,6,n);
%             plot(time,zscore(meg(iir,:)),'Linewidth',2)
%             hold on
%             %             plot(time,X.Estimate,'k')
%             %             plot(time(ind),X.Estimate(ind), '*r')
%             
%             plot(time,x.tStat,'k')
%             plot(time(ind),x.tStat(ind), '.r')
% 
%             if nrois == length(aal_labels)
%                 title(aal_labels{iir})
%             elseif nrois == length(channels)
%                 title(channels{iir})
%             end
% %             xlabel('time (s)'); ylabel('tStat')
%             grid on;
% %             sb.XTick = -0.2:0.1:1;
% %             sb.XTickLabel([1,2,4,5,7,8,10,11,13]) = {''};
% %             sb.YTick = -6:2:6;
%             xlim(timew); ylim([-5 5])
% %             sb.XTickLabel = sb.XTick;
%         end
%             
%     end
% %     if n ~= 0
% %         saveas(gcf,sprintf('~/matlab/figures/%s_%s.tif',freq,fit_parameters{m}))
% %     end
% end


%% Calculate Threshold free cluster enhancement
close all
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
clusternull = cell(5,1);
condition = freq((ind+1):end);
for m = 1:5
    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s_%s/lmixmodel_%s',data_path,freqb,condition,fit_parameter))
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
% snull = sort(abs(clusternull(:)));
% Nnull = size(snull,1);
% uplim =  snull(ceil(Nnull*(1-alpha)),1);
% lowlim = -uplim;
% 
% 
snull = sort(clusternull);
Nnull = size(snull,1);
lowlim = snull(floor(Nnull*alpha),2);
uplim =  snull(ceil(Nnull*(1-alpha)),1);

dx = 1;
for m = 1:length(fit_parameters)
    
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
%         saveas(gcf,sprintf('~/matlab/figures/%s_%s.tif',freq,fit_parameters{m}))
    end

end
return
%%
fit_parameter = 'E_LTA';
meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
meg = reshape(meg, npoints,nrois,size(meg,2));

m = find(strcmp(fit_parameters,fit_parameter));

opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);
Xr = eval(['X.',fit_parameter]);

t = time > -1 & time < -0.4;
iir = find(strcmp(aal_labels,'FrontalSupOrb_L'));

megr = squeeze(meg(t,iir,:));
figure
scatter(Xr,mean(megr,1));
xlabel(fit_parameter); ylabel(freq)
%%
% 
% %%
% title_parameters = fit_parameters;
% title_parameters{strcmp(fit_parameters,'LTA_sum')} = '\SigmaLTA';
% title_parameters{strcmp(fit_parameters,'RPE_sum')} = '\SigmaRPE';
% % close all
% sourcemni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');
% % sourcemni = ft_read_mri('~/MNI152_T1_2009c.nii');
% 
% atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
% sourcemni.tfce = NaN(size(sourcemni.anatomy));
% m =1;
% 
% % centroid of ROIs
% locs = zeros(length(param_list),3);
% for iir  = 1:length(param_list)
%     ind = find(atlas.tissue == iir);
%     [x,y,z] = ind2sub(size(atlas.tissue),ind);
%     locs(iir,:)  = mean([x,y,z]); % centroid        
% end
% % Find relative distance of ROI
% D = zeros(length(param_list));
% for iir  = 1:length(param_list)
%     % distance in voxels
%     D(iir,:) = sqrt(sum((locs(iir,:) - locs).^2,2));
% end
% 
% for t = 0.9 %-.15:0.05:0.95
% twind = t+0.005*[-1 1];
% tpoints = time>= twind(1) & time <= twind(2);
% 
% for iir = 1:length(param_list)
%     % keep TFCE value
% %     sourcemni.tfce(atlas.tissue==iir) = abs(mean(clusteraal{m}(iir,tpoints),2));
%     % Binary yes/no 
%     sourcemni.tfce(atlas.tissue==iir) = ...
%         any(clusteraal{m}(iir,tpoints) > uplim) | ...
%         any(clusteraal{m}(iir,tpoints) < lowlim)  ;
% end
% 
% % crang = [uplim 60];
% crang = [0 2.5];
% cfg = [];
% cfg.method        = 'slice';
% % cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
% cfg.funparameter = 'tfce';
% cfg.maskparameter = 'tfce';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourcemni);
% title(sprintf('%s : time %.0f : %.0f ms',title_parameters{m},twind(1)*1e3,twind(2)*1e3))
% % saveas(gcf,sprintf('~/matlab/figures/%s_%s%.0f.tif',freq,fit_parameters{m},twind(1)*1000))
% 
% end
% 

