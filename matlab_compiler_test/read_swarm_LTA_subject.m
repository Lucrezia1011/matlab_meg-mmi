clear all 
close all
clc
%% Read Swarm output

% stype = 'sensors';
stype = 'aal_mu5max_Z';

if strcmp(stype,'sensors')
    nrois = 269;
else
    nrois = 116;
end
freq  = ['evoked_',stype];


data_path = '/data/MBDU/MEG_MMI3/results/mmiSub_level_LTA/';

% nrois = 269;
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/';


latent_vars_name = 'latent_vars.csv';
opts = detectImportOptions([data_path,latent_vars_name]);
X = readtable([data_path,latent_vars_name],opts);

fit_parameters = X.Properties.VariableNames(2:(end-3));


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


dataty = 'posneg';

if strcmp(dataty,'pos')
    dd = 1;
    npoints = 1200;
    meg_data_name = ['positive_',stype,'.txt'];
elseif strcmp(dataty,'neg')
    dd = 2;
    npoints = 1200;
    meg_data_name = ['negative_',stype,'.txt'];
elseif strcmp(dataty,'posneg')
    dd = 3;
    npoints = 1200;
    meg_data_name = ['posneg_',stype,'.txt'];
end
        
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


if ~exist('meg','var')
    meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
    meg = mean(meg,2);
    meg = reshape(meg, nrois,npoints);
end
time = linspace(-3,3,npoints);

%% Check for missing points
% 
meg_data_names = {['positive_',stype,'.txt'];['negative_',stype,'.txt'];...
    ['posneg_',stype,'.txt']};


n = 0;
runcompiled = 'run_mmi_LTA_subjects.sh';               

Nnames = {'Npos';'Nneg';'Nposneg'};
outname = {[freq,'_pos'];[freq,'_neg'];[freq,'_posneg']};
% 
% 
% command_list = [];
% for d = 1:length(meg_data_names)
% for m = 1:length(fit_parameters)
%     fit_parameter = fit_parameters{m};
%     cd(sprintf('%s/%s/glmmodel_%s',data_path,outname{d},fit_parameter))
%     for ii = 1:length(param_list)
%         filename = sprintf('ROI_%s.csv',param_list{ii});
%         
%         if ~exist(filename,'file')
%             n = n +1;
%             command_list = [command_list ...
%                 sprintf(['export MCR_CACHE_ROOT=/lscratch/$SLURM_JOB_ID; '...
%                 '  test -d /lscratch/$SLURM_JOB_ID/v96 || tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v96.tar.gz '...
%                 '  && ~/matlab/matlab_compiler_test/%s '...
%                 ' /lscratch/$SLURM_JOB_ID/v96 %s %s %s %s %s %s %s %s\n'],runcompiled,...
%                 meg_data_names{d},latent_vars_name,param_list{ii},num2str(npoints),Nnames{d},fit_parameter,outname{d},data_path)];
%             
%         end
%         
%     end
% end
% end


%
% % identify missing parameters: run as normal function
% cd ~/matlab/matlab_compiler_test/
% fit_parameter = 'b';
% param_list = {'083'};
% d = 1;
% parfor ii = 1:n
%     mmi_LTA_subjects(meg_data_names{d},latent_vars_name,param_list{ii},num2str(npoints),Nnames{d},fit_parameter,outname{d},data_path)
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
for m = [2:5,9:length(fit_parameters)]
    X = [];

    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s/glmmodel_%s',data_path,outname{dd},fit_parameter))
    
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

%% plots
timex = linspace(-.5,1,200*1.5);

for m = [2:5,9:length(fit_parameters)]
    fit_parameter = fit_parameters{m};
    X = Xfit{m};
    [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
    n = 0;
    for ii = 1:length(param_list)
        x = X(r==ii & t>=500 & t<800,:);
%         x = X(r==ii,:);
        ind = find(x.pValue<(0.05));
        if any(ind)
            if n == 0
                figure(m); clf
                set(gcf,'color','w','name',fit_parameter,'position',[637   204   859   661])
            end
            n = n+1;
            subplot(4,4,n)
            plot(time,zscore(meg(ii,:)),'Linewidth',2)
            hold on
            %             plot(time,X.Estimate,'k')
            %             plot(time(ind),X.Estimate(ind), '*r')
            
            plot(timex,x.tStat,'k')
            plot(timex(ind),x.tStat(ind), '*r')
            if nrois == length(aal_labels)
                title(aal_labels{ii})
            elseif nrois == length(channels)
                title(channels{ii})
            end
            xlabel('time (s)'); ylabel('tStat')
             xlim([-.5 1])
        end
    end
    
end

%% cluster test

npointsx = 200*1.5;
timex = linspace(-.5,1,npointsx);

clusteraal = cell(1,length(fit_parameters));
dx = 1;
hexp = 2; % 2 by default
Eexp = 0.5;
for m =  [2:5,9:length(fit_parameters)]
    fit_parameter = fit_parameters{m};
    X = Xfit{m};
    [r,tt] = ind2sub(size(meg),X.index); % find real ROI positions
    n = 0;
    clusteraal{m} = zeros(length(param_list),npointsx);
    for iir = 1:length(param_list)
        
        TFCE = zeros(1,npointsx);

        x = X(r==iir & tt>=500 & tt<800 ,:);
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
            a{end} = ii(end)+1:npointsx;
        else
            a = cell(1,length(ii)/2+1);
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)-1
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
            a{end} = ii(end)+1:npointsx;
        end
        
        
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
                    eE = (sum(e)*dx)^Eexp;
                    
                    A = A + eE*(hp^hexp)*dh;
                    hp = hp + dh;
                end
                TFCE(p+a{ii}(1)-1) = A;
            end
            
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
            a{end} = ii(end)+1:npointsx;
        else
            a = cell(1,length(ii)/2+1);
            a{1} = 1:ii(1);
            for jj = 1:floor(length(ii)/2)-1
                a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
            end
            a{end} = ii(end)+1:npointsx;
        end
        
        
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
                    eE = (sum(e)*dx)^Eexp;
                    
                    A = A + eE*(hp^hexp)*dh;
                    hp = hp + dh;
                end
                TFCE(p+a{ii}(1)-1) = -A;
            end
            
        end
        
        ind0 = find(x.pValue<(0.05/npointsx));
        ind = find(abs(TFCE)>100);
        if any(ind)  
            
            if n == 0
                figure(m); clf
%                 set(gcf,'color','w','name',fit_parameter,'position',[637   204   859   661])
%                 set(gcf,'color','w','name',fit_parameter,'position',[345   335   1301   510])
                set(gcf,'color','w','name',fit_parameter,'position',[ 345  215   1094 630])
            end
            n = n+1;
            subplot(4,5,n)
            plot(time,zscore(meg(iir,:)),'Linewidth',2)
            hold on
            %             plot(time,X.Estimate,'k')
            %             plot(time(ind),X.Estimate(ind), '*r')
            
            plot(timex,x.tStat,'k')
            plot(timex(ind),x.tStat(ind), '.r')
            xlim([-.5 1])
            if nrois == length(aal_labels)
                title(aal_labels{iir})
            elseif nrois == length(channels)
                title(channels{iir})
            end
            xlabel('time (s)'); ylabel('tStat')
        end
            
        clusteraal{m}(iir,:) = TFCE;
    end
end


%%
title_parameters = fit_parameters;

close all
sourcemni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii');
atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
sourcemni.tfce = NaN(size(sourcemni.anatomy));
m =3;

for t = -0.2:0.05:0.7
twind = t+[0.0 0.1];
tpoints = time>= twind(1) & time <= twind(2);

for iir = 1:length(param_list)
    sourcemni.tfce(atlas.tissue==iir) = abs(mean(clusteraal{m}(iir,tpoints),2));
end

crang = [5 35];
cfg = [];
cfg.method        = 'slice';
% cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
cfg.funparameter = 'tfce';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourcemni);
title(sprintf('%s : time %.2fs-%.2fs',title_parameters{m},twind(1),twind(2)))
% saveas(gcf,sprintf('~/matlab/figures/%s_%s%d.tif',freq,fit_parameters{m},twind(1)*100))

end



