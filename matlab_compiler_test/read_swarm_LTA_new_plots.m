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

if strcmp(freq,'evoked_outcome')
    npoints = 360;
    timew = [-.2,1];
elseif strcmp(freq,'evoked_choice')
    npoints = 390;
    timew = [-.5,.8];
else
    npoints = 360;
    timew = [-3,3];
end
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
Xv = readtable([data_path,latent_vars_name],opts);

fit_parameters = Xv.Properties.VariableNames(3:end);
% fit_parameters([6,7])= [];


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

time = linspace(timew(1),timew(2),npoints);

subs = unique(Xv.subject);

    meg = dlmread(sprintf('%s%s',data_path,meg_data_name));
    meg = reshape(meg, npoints,nrois,size(meg,2));
   
 
%     timer = time>=0.25 & time <= .35;
%     
%     for ii = 110; clf;%116; clf
%         
        
%     nneg = Xv.RPE<-3; % 21
%     npos = Xv.RPE<0.726 & Xv.RPE>-0.726;
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
%     ylim([-0.7 0.7]); grid on; legend('RPE=0','RPE<-3')
%     pause(1)
%     
%     
%     subplot(212)
%     npos = Xv.RPE>5;
%     nneg = Xv.RPE<0.87 & Xv.RPE>-0.88;
%      plot(time, mean(meg(:,ii,npos),3))
%     hold on
%     plot(time, mean(meg(:,ii,nneg),3)); title(aal_labels{ii})
%     fill([time,fliplr(time)],[mean(meg(:,ii,npos),3)+std(meg(:,ii,npos),0,3)/sqrt(nnz(npos));...
%         flipud(mean(meg(:,ii,npos),3)-std(meg(:,ii,npos),0,3)/sqrt(nnz(npos)))],...
%         [0 0 1],'facealpha',0.3,'edgecolor','none')
%     fill([time,fliplr(time)],[mean(meg(:,ii,nneg),3)+std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg));...
%         flipud(mean(meg(:,ii,nneg),3)-std(meg(:,ii,nneg),0,3)/sqrt(nnz(nneg)))],...
%         [1 0 0],'facealpha',0.3,'edgecolor','none');
%     ylim([-0.7 0.7]); grid on; legend('RPE>5','RPE=0')
    
%     end

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

%%
mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
    % mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii');
    ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
    gridres = 5;
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
    sourcemodel.coordsys = 'mni';
    atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
    
    atlas = ft_convert_units(atlas,sourcemodel.unit);
    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = 'tissue';
    sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

%%    
close all
for tm = 0.05:0.05:0.65
    tt = time>=(tm-0.05) & time <= (tm+0.05);
%     tt = time>=-0.1 & time <= 0;  
%     tt = time>=0.1 & time <= 0.2;
    
    Tstat  = zeros(size(sourcemodel.inside));
    for ii = 1:116
        ind = find(sourcemodelAAL.tissue == ii);
        Tstat(ind) = mean(abs(meg(ii,tt)),2);
    end
    sourceant =[];
    sourceant.pow = Tstat;
    sourceant.dim = sourcemodel.dim;
    sourceant.inside = sourcemodel.inside;
    sourceant.pos = sourcemodel.pos;
    cfg = [];
    cfg.parameter = 'pow';
    sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
    
    crang = [0.08 0.3];
    cfg = [];
    cfg.method        = 'slice'; %'ortho'
    cfg.funparameter = 'pow';
    cfg.maskparameter = 'pow';
    cfg.funcolormap  = 'auto';
    cfg.funcolorlim   = crang;
    cfg.opacitylim = crang;
    
    ft_sourceplot(cfg, sourceout_Int);
    title(sprintf('Evoked response %.0f-%.0fms',(tm-0.05)*1000,(tm+0.05)*1000))
end