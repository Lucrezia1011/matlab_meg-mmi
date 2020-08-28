clear all
close all
% clc
meginfo = readtable('~/MEG_participantsinfo.csv');
% make a list of excluded recordings (too short <5min of data)
data_exclude = {'sub-24201_task-mmi3_run-1_meg.ds';...
    'sub-22694_task-mmi3_run-2_meg.ds'; ...
    'sub-22694_task-mmi3_run-3_meg.ds'; ...
    'sub-23999_task-mmi3_run-3_meg.ds';...
    'sub-22812_task-mmi3_run-2_meg.ds';...
    'sub-22658_task-mmi3_run-1_meg.ds'};

data_list = [];


% exclude subject 10: did not perform task correctly
% subject 24 : metal artefacts
% subjects 26,49,53: no co-registration
Nlist = 1:56;
Nlist([10,24]) = []; %Nlist([10,24,26,49,53]) = [];
zz= 0;
for sn = Nlist %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sdan = num2str(meginfo.SDAN(sn));
    cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            zz = zz +1;
            data_list{zz} = data_name;
        end
    end
    
   
end

%%
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca
% 
% % roiopt = 'c';
% roiopt = 'grid';
% % tfsopt = 'pwelch';
% tfsopt = 'm';
% 
gridres= 5;

% freql=[ 4 8; 8 13; 13 25; 25 40; 40 150]; filter_type = 'but';
% freql = [1 4]; filter_type = 'firls';

roiopt = 'grid'; % running for grid
freqband = [25 40]; 
for ii = 1:length(data_list)
    data_name = data_list{ii};
    sub = data_name(5:9);
    mmiPreMoodPower(data_name,roiopt,gridres,freqband); % pre mood
    
end

%
% 
% if strcmp(roiopt,'grid') || strcmp(roiopt,'sens')
%     for freq = 1:size(freql,1)
%         for ii = 1:length(data_list)       
%             mmiPreMoodPower(data_list{ii},roiopt,gridres,freql(freq,:)); % pre mood
% %             mmi_grid_prep_PowerITI(param_list{ii},roiopt,gridres,freql(freq,:))
% %             mmi_grid_prep_Power(param_list{ii},roiopt,gridres,freql(freq,:)) % pre mood
% %             mmi_grid_prep_PowerA(param_list{ii},roiopt,gridres,freql(freq,:),filter_type)
% %             % anticipation period          
%         end
%     end    
% else
%     for ii = 1:length(data_list)
%         mmi_grid_prep(data_list{ii},roiopt,gridres,tfsopt)
%     end
% end

% resave ltv variables from delta

return

%% Find all common gradiometers

BadChannels = [];
for s = 1:length(data_list)
    
    data_name = data_list{s};
    sub = data_name(5:9);
    cd(['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'])
    
    % Get Bad channel names
    fid = fopen([data_name,'/BadChannels']);
    BadChannel = textscan(fid,'%s');
    fclose(fid);
    BadChannels = cat(1,BadChannels,BadChannel{1});
    
end

BadChannels = unique(BadChannels);
% get MEG channel names
hdr = ft_read_header(data_name);
channels = hdr.label(strcmp(hdr.chantype,'meggrad'));
% Delete Bad channels
chanInd = zeros(size(channels));
for iiC = 1:length(BadChannels)
    chanInd = chanInd | strcmp(channels,BadChannels{iiC});
end
channels(chanInd) = [];

save('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors','channels');
%%
Pall = [];
Vall = [];

ltvMood = [];

Ytfs = [];
Ytfsp = [];

if strcmp(roiopt,'sens')

    if ~exist('sensall','var')
        load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    end

    cd(['/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/'])
    hdr = ft_read_header(data_list{1});
    channelsall = hdr.label(strcmp(hdr.chantype,'meggrad'));
end

for sn = 1:length(data_list) 
    
    
    data_name = data_list{sn};
    sub = data_name(5:9);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
    
    
    if strcmp(roiopt,'grid')
        cd(data_path)
        load('leadfields_5mm.mat')
        if exist('gridall','var')
            gridall = gridall & grid.inside;
        else
            gridall = grid.inside;
        end
        
        load(sprintf('pre_mood_grid_%.0f-%.0fHz',...
            freqband(1),freqband(2)))
        
        Pmni = zeros(size(grid.inside,1),size(P,2));
        Pmni(grid.inside,:) = P;
        Pall = cat(2,Pall,Pmni);
        
        
    elseif  strcmp(roiopt,'sens')
        
        cd(['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'])
        
        % Get Bad channel names
        fid = fopen([data_name,'/BadChannels']);
        BadChannel = textscan(fid,'%s');
        BadChannel = BadChannel{1};
        fclose(fid);
        channelSub = channelsall;
        % Delete Bad channels
        chanInd = zeros(size(channelsall));
        for iiC = 1:length(BadChannel)
            chanInd = chanInd | strcmp(channelsall,BadChannel{iiC});
        end
        channelSub(find(chanInd)) = [];
        cd(data_path)
        
        load(sprintf('pre_mood_sens_%.0f-%.0fHz',...
            freqband(1),freqband(2)))
        
        [~,~,ind]= intersect(channels,channelSub);
        Vall = cat(2,Vall,V(ind,:));
        
    else
        load(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/',sub,'_',num2str(n)])
        nrois = length(TFS);
        
        [nf,nt] = size(TFS{1});
        TFS = cell2mat(TFS);
        TFS = reshape(TFS,[nf,nt,nrois]);
        Ytfs = cat(2,Ytfs,TFS);
        
        [nf,nt] = size(TFSp{1});
        TFSp = cell2mat(TFSp);
        TFSp = reshape(TFSp,[nf,nt,nrois]);
        Ytfsp = cat(2,Ytfsp,TFSp);
    end
    
    ltvmood.recording = repmat(sn,size(ltvmood,1),1);
    if isempty(ltvMood)
        ltvMood = ltvmood;
    else
        ltvMood(end+(1:size(ltvmood,1)),:) = ltvmood;
    end
    
    
end

clear ltvmood TFS P Pmni grid
if strcmp(roiopt,'grid')
    Pall = Pall(gridall,:); 
end
%% Write data (grid)
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';
dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
    freqband(1),freqband(2)),Pall);

dlmwrite([out_path,'/mni_grid.txt'],gridall);
writetable(ltvMood,[out_path,'/latent_vars.csv']);
  
%% Write data (sens)
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood';
dlmwrite(sprintf('%s/powersens_%.0f-%.0fHz.txt',out_path,...
    freqband(1),freqband(2)),Vall);
writetable(ltvMood,[out_path,'/latent_vars.csv']);


%% Plot Linear mixed effects model for grid
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';
ltvMood = readtable([out_path,'/latent_vars.csv']);

freql=[ 25 40];


mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';


datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/powergrid_%.0f-%.0fHz',...
    freql(1),freql(2));

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/mni_grid.txt');
nf = size(freql,1);
Tmood = zeros(size(gridall));
Te = zeros(size(gridall));
Tesum = zeros(size(gridall));

param_list = cell(1,15);
for nn = 1:14
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

cd([datapath,'/lme_mood'])

X = zeros(nnz(gridall),1);
for nn = 1:14
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
end
Tmood(gridall==1) = X;


cd([datapath,'/lme_E_sum'])

X = zeros(nnz(gridall),1);
for nn = 1:14
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
end

Tesum(gridall==1) = X;


cd([datapath,'/lme_E'])

X = zeros(nnz(gridall),1);
pV =  zeros(nnz(gridall),1);
for nn = 1:14
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
end
Te(gridall==1) = X;

fprintf(sprintf('Loaded frequency %.0f-%.0fHz\n',freql(1),freql(2)))



%% Plot p-value
figure; set(gcf,'color','w')

p = sort(pV);
semilogy(p)
hold on
semilogy(0.05./(nnz(gridall):-1:1))
%     grid on
xlabel('voxels')
ylabel('p-value')
legend('p-values','FDR','location','best')
N = nnz(p'<0.05./(nnz(gridall):-1:1));

title(sprintf('E: p-value of %.0f voxels < 0.05 (FDR)',N))
%%

for p =  2:3
    sourceant =[];

    switch p
        case 1
            sourceant.pow = Tmood;
            fit_parameter = 'mood';
        case 2
            sourceant.pow = Te;
            fit_parameter = 'E';
        case 3
            sourceant.pow = Tesum;
            fit_parameter = 'E sum';
    end
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';


crang = [2 5];
cfg = [];
cfg.method        = 'ortho'; %'ortho'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
title(sprintf('%s : %.0f-%.0fHz\npeak t-value %.1f',...
    fit_parameter,freql(1),freql(2),max(abs(sourceant.pow(:)))))

saveas(gcf,sprintf('~/matlab/figures/%s_gamma_slice.png',fit_parameter))
% saveas(gcf,sprintf('~/matlab/figures/ELTAgamma_peak1.png'))
end

%% Plot source data
% megP = dlmread([datapath,'.txt']);

% for ii = 1:max(ltvMood.recording)
%     M = megP(:,ltvMood.recording==ii);
%     megP(:,ltvMood.recording==ii) = (M) / std(M(:));
% end
% meg = zeros(size(sourcemodel.inside));
% meg(gridall==1) = mean(megP,2); % simple mean


ii = 6;

meg = zeros(size(sourcemodel.inside));
M = megP(:,ltvMood.recording==ii);
meg(gridall==1) = mean(M,2); 

sourceant.pow = meg;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';


crang = [];
cfg = [];
cfg.method        = 'slice'; %'ortho'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);

%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';
freql=[1 4; 4 8; 8 13; 13 25; 25 40; 40 150];
if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
end

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/';
ff = 5;


freq = sprintf('powersens_%.0f-%.0fHz',...
    freql(ff,1),freql(ff,2));

meg_data_name = sprintf('%s.txt',freq);
meg = dlmread([data_path,meg_data_name]);

freq = sprintf('powersens_%.0f-%.0fHz_rec-effect',...
    freql(ff,1),freql(ff,2));
outpath = sprintf('%s%s/',data_path,freq);
nn =1;
fit_parameters = {'mood';'E';'E_sum';'RPE';'RPE_sum'};
Tfit = cell(1,length(fit_parameters));
pfit = cell(1,length(fit_parameters));
for ii = 1:length(fit_parameters)
    cd([outpath,'lme_',fit_parameters{ii}])
    opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
    Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
    Tfit{ii} = Xv.tStat;
    pfit{ii} = Xv.pValue;
end


T = struct;
T.label = channels;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position',[176 748 1262 800])

for ii = 1:length(fit_parameters)
subplot(2,3,ii)
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
subplot(236)
caxis(cfg.zlim)
colorbar

% saveas(gcf,sprintf('~/matlab/figures/pre_mood_%s.png',freq))

% plot the MEG data
ltvMood = readtable('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/latent_vars.csv');
subs = unique(ltvMood.subject);
megz = zeros(size(meg,1),length(subs));
for ss = 1:length(subs)
    megm = mean(meg(:,ltvMood.subject==subs(ss)),2); % mean over trials
    megz(:,ss)  = zscore(megm);
end
figure(ff+2); clf; set(gcf,'color','w')
T.trial{1} = mean(megz,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = channels;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-1 1];
cfg.colorbar = 'yes';
ft_topoplotER(cfg, T);
title(sprintf('Grand-average (Z-scored) of 25-40Hz variance\nduring mood rating preparation'))
% saveas(gcf,sprintf('~/matlab/figures/pre_mood_%s_avg.png',freq))

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


%% LME for TFS analysis
lme_formula2 = 'MEG ~ E + (E|subject) + (1|trial) + (1|recording)';
lme_formula1 = 'MEG ~ mood +(mood|subject) + (1|trial) + (1|recording)';
lme_formula3 = 'MEG ~ RPE + (RPE|subject) + (1|trial) + (1|recording)';

if strcmp(tfsopt,'m')
% Run linear mixed model with frequency power envelope
% Multitaper
    fskip = find(isnan(Ytfs(:,1,1)));
    Ytfs(fskip,:,:) = [];
    freqs(fskip) = [];
    nf = length(freqs);
elseif strcmp(tfsopt,'p')
% pwelch 
    freql=[1 4; 4 8; 8 13; 13 25; 25 40; 40 150];
    nf = size(freql,1);
    ytfs = zeros(nf,size(Ytfsp,2),nrois);

    for ff = 1:size(freql,1)
        ytfs(ff,:,:) = mean(Ytfsp(F>=freql(ff,1) & F<=freql(ff,2),:,:),1);
    %     ytfs(ff,:,:) = max(Ytfsp(F>=freql(ff,1) & F<=freql(ff,2),:,:),[],1);
    end

end

tfs_coeff1 = cell(nf,nrois);
tfs_coeff2 = cell(nf,nrois);
tfs_coeff3 = cell(nf,nrois);


% what if I tried to find a fixed effect of trial?
parfor n = 1:nrois
    X = ltvMood;
    
    for ff = 1:nf
%         X.MEG = Ytfs(ff,:,n)';
        X.MEG = ytfs(ff,:,n)';
        lme = fitlme(X,lme_formula1); 
        tfs_coeff1{ff,n} = [lme.Coefficients.tStat(2); ...
            lme.Coefficients.pValue(2)]; 
        lme = fitlme(X,lme_formula2); 
        tfs_coeff2{ff,n} = [lme.Coefficients.tStat(2); ...
            lme.Coefficients.pValue(2)]; 
        lme = fitlme(X,lme_formula3); 
        tfs_coeff3{ff,n} = [lme.Coefficients.tStat(2); ...
            lme.Coefficients.pValue(2)]; 
    end
    clc; fprintf('Done roi %.0f/%.0f\n',n,nrois);
end

tfs_coeff1 = cell2mat(tfs_coeff1);
tfs_coeff1 = reshape(tfs_coeff1,[2,nf,nrois]);

tfs_coeff2 = cell2mat(tfs_coeff2);
tfs_coeff2 = reshape(tfs_coeff2,[2,nf,nrois]);

tfs_coeff3 = cell2mat(tfs_coeff3);
tfs_coeff3 = reshape(tfs_coeff3,[2,nf,nrois]);

%%
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');
freqnames = {'\delta';'\theta';'\alpha';'\beta';'Low-\gamma';'High-\gamma'};


p = 0.01;

% figure; pcolor(freqs,1:nrois,squeeze(tfs_coeff(1,:,:))')
% shading interp; colorbar; caxis([-4,4])


figure; set(gcf,'color','w'); imagesc(squeeze(tfs_coeff1(1,:,:))')
colorbar; caxis([-4,4]); colormap jet
set(gca,'XTick',1:nf,'XTickLabel',freqnames); ylabel('ROI')
title('Mood')
[ff,nn] = find(squeeze(tfs_coeff1(2,:,:))<p );

clc
for ii = 1:length(nn)
    fprintf('Mood: %s, %s, T-stat %.2f\n',...
       aal_labels{nn(ii)},freqnames{ff(ii)}, tfs_coeff1(1,ff(ii),nn(ii)))
end


figure;  set(gcf,'color','w'); imagesc(squeeze(tfs_coeff2(1,:,:))')
colorbar; caxis([-4,4]); colormap jet
set(gca,'XTick',1:nf,'XTickLabel',freqnames); ylabel('ROI')
title('LTA Expectation');
[ff,nn] = find(squeeze(tfs_coeff2(2,:,:))<p );


for ii = 1:length(nn)
    fprintf('E:  %s, %s, T-stat %.2f\n',...
        aal_labels{nn(ii)},freqnames{ff(ii)}, tfs_coeff2(1,ff(ii),nn(ii)))
end

figure;  set(gcf,'color','w');imagesc(squeeze(tfs_coeff3(1,:,:))')
colorbar; caxis([-4,4]); colormap jet
set(gca,'XTick',1:nf,'XTickLabel',freqnames); ylabel('ROI')
title('LTA RPE');
[ff,nn] = find(squeeze(tfs_coeff3(2,:,:))<p );


for ii = 1:length(nn)
    fprintf('RPE: %s, %s, T-stat %.2f\n',...
        aal_labels{nn(ii)},freqnames{ff(ii)}, tfs_coeff3(1,ff(ii),nn(ii)))
end
%%

pv = 0.05;
M = nrois;
tlim = 6;

C = corr(YDelta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

figure; set(gcf,'color','w','name',lme_formula) 

subplot(2,3,1); 
semilogy(tfs_coeff(2,:),'k'); 
hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(tfs_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
xlabel('ROI');
title('Delta power'); grid on

aal_labels(tfs_coeff(2,:)<=alpha)


C = corr(YTheta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);


subplot(2,3,2); 
semilogy(theta_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(theta_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Theta power'); grid on
xlabel('ROI');

aal_labels(theta_coeff(2,:)<=alpha)

C = corr(YAlpha');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,3); 
semilogy(alpha_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r')
ylabel('p-value')
yyaxis right; plot(alpha_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Alpha power'); grid on
xlabel('ROI');

aal_labels(alpha_coeff(2,:)<=alpha)

C = corr(YBeta');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,4);
semilogy(beta_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r');
ylabel('p-value')
yyaxis right; plot(beta_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Beta power'); grid on
xlabel('ROI');

aal_labels(beta_coeff(2,:)<=alpha)

C = corr(YGamma');
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

subplot(2,3,5); 
semilogy(gamma_coeff(2,:),'k'); hold on
semilogy([0,nrois+1],[1 1]*alpha,'r');
ylabel('p-value')
yyaxis right; plot(gamma_coeff(1,:))
ylabel('T-stat'); ylim(tlim*[-1 1])
title('Gamma power'); grid on
xlabel('ROI');

aal_labels(gamma_coeff(2,:)<=alpha)