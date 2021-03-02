clear all
close all
clc

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn = 1:16 %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 1.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') ...
                && exist([data_path,name_list(ii).name,'/beamforming/ICA_artifacts.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)   
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

%%
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')
addpath ~/ppyll1/matlab/svdandpca


% roiopt = 'c';
roiopt = 'grid';
% tfsopt = 'pwelch';
tfsopt = 'm';

gridres= 5;

% freql=[ 4 8; 8 13; 13 25; 25 40; 40 150]; filter_type = 'but';
% freql = [1 4]; filter_type = 'firls';

roiopt = 'sens';
freql = [25 40]; filter_type = 'but';


if strcmp(roiopt,'grid') || strcmp(roiopt,'sens')
    for freq = 1:size(freql,1)
        for ii = 1:length(param_list)
            mmi_grid_prep_Powerltv(param_list{ii},roiopt,gridres,freql(freq,:)); % pre mood
%             mmi_grid_prep_PowerITI(param_list{ii},roiopt,gridres,freql(freq,:))
%             mmi_grid_prep_Power(param_list{ii},roiopt,gridres,freql(freq,:)) % pre mood
%             mmi_grid_prep_PowerA(param_list{ii},roiopt,gridres,freql(freq,:),filter_type)
%             % anticipation period          
        end
    end    
else
    for ii = 1:length(param_list)
        mmi_grid_prep(param_list{ii},roiopt,gridres,tfsopt)
    end
end

% resave ltv variables from delta

return

%% Find all common gradiometers
outpath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback/';

for s = 1:length(param_list)
    data_name = param_list{s};
    sub = data_name(1:5);
    
    sens = ft_read_sens(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name],...
        'senstype','meg');
    if s == 1
        sensall = sens.label(strcmp(sens.chantype,'meggrad'));
    else
        sensall = intersect(sensall,sens.label(strcmp(sens.chantype,'meggrad')));
    end
end

save('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors','sensall');
%%
% out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/ITI';
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood';
% out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_feedback';

Pall = cell(1,size(freql,1));
Vall = cell(1,size(freql,1));

ltvMood = [];

Ytfs = [];
Ytfsp = [];
z = 0;
for sn = 1:16
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') ...
                && exist([data_path,name_list(ii).name,'/beamforming/ICA_artifacts.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)  
        z = z+1;
        data_name = data_names{runs};
        n = str2double(data_name(end-3));
        if isnan(n)
            n = 1;
        end
        
        
        if strcmp(roiopt,'grid')
            load([data_name,'/beamforming/leadfields_5mm.mat'])
            if exist('gridall','var')
                gridall = gridall & grid.inside;
            else
                gridall = grid.inside;
            end
            
            for freq = 1:size(freql,1)
                load(sprintf('%s/%.0f-%.0fHz_%s_%.0f',...
                    out_path,freql(freq,1),freql(freq,2),sub,n))
                Pmni = zeros(size(grid.inside,1),size(P,2));
                Pmni(grid.inside,:) = P;
                Pall{freq} = cat(2,Pall{freq},Pmni);
            end
        
        elseif  strcmp(roiopt,'sens')
            if ~exist('sensall','var')
                load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat')
            end
            for freq = 1:size(freql,1)
            load(sprintf('%s/%.0f-%.0fHz_%s_%.0f',...
                    out_path,freql(freq,1),freql(freq,2),sub,n))
            sens = ft_read_sens(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name],...
                'senstype','meg');
            [~,~,ind]= intersect(sensall,sens.label(strcmp(sens.chantype,'meggrad')));
            Vall{freq} = cat(2,Vall{freq},V(ind,:));
            end
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
        
        ltvmood.recording = repmat(z,size(ltvmood,1),1);
        if isempty(ltvMood)
            ltvMood = ltvmood;  
        else
            ltvMood(end+(1:size(ltvmood,1)),:) = ltvmood;
        end
                
    end
end

clear ltvmood TFS P Pmni grid
for freq = 1:size(freql,1)
    Pall{freq} = Pall{freq}(gridall,:);
end
%% Write data (grid)
for freq = 1:size(freql,1)
    dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
                    freql(freq,1),freql(freq,2)),Pall{freq});
end
dlmwrite([out_path,'/mni_grid.txt'],gridall);
writetable(ltvMood,[out_path,'/latent_vars.csv']);
  
%% Write data (sens)
for freq = 1:size(freql,1)
    dlmwrite(sprintf('%s/powersens_%.0f-%.0fHz.txt',out_path,...
                    freql(freq,1),freql(freq,2)),Vall{freq});
end
writetable(ltvMood,[out_path,'/latent_vars.csv']);

%% Read data
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood';
ltvMood = readtable([out_path,'/latent_vars.csv']);
freql=[1 4; 4 8; 8 13; 13 25; 25 40; 40 150];
for freq = 5%1:5
% P = dlmread(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
%                     freql(freq,1),freql(freq,2)));
              
subs = unique(ltvMood.subject);
Pdiff = zeros(size(P,1),length(subs));
for s=1:length(subs)
    Ps = P(:,ltvMood.subject==subs(s));
    Esum = ltvMood.E(ltvMood.subject==subs(s));
    Pdiff(:,s) = (mean(Ps(:,Esum>median(Esum)+2),2) - mean(Ps(:,Esum<median(Esum)-2),2))./mean(Ps,2);

end

Pall = zeros(size(gridall));
% Pall(gridall==1) = mean(P(:,ltvMood.recording==2),2);
Pall(gridall==1) = mean(Pdiff(:,1:end),2);

sourceant.pow = Pall;
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
title(sprintf('%.0f-%.0fHz',freql(freq,1),freql(freq,2)))
end
%% Plot Linear mixed effects model for grid
ltvMood = readtable([out_path,'/latent_vars.csv']);

freql=[1 4; 4 8; 8 13; 13 25; 25 40; 40 150];


mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/mni_grid.txt');
nf = size(freql,1);
Tmood(1:nf) = {zeros(size(gridall))};
Te(1:nf) = {zeros(size(gridall))};
Trpe(1:nf) = {zeros(size(gridall))};

param_list = cell(1,15);
for nn = 1:15
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

for freq = 5%1:nf
    datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/powergrid_%.0f-%.0fHz/',...
        freql(freq,1),freql(freq,2));
    cd([datapath,'lme_mood'])
    
    X = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    end
    Tmood{freq}(gridall==1) = X;
    
    
    cd([datapath,'lme_Esum'])
    
    X = zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
    end
    
    Trpe{freq}(gridall==1) = X;
    
    
    cd([datapath,'lme_E'])
    
    X = zeros(nnz(gridall),1);
    pV =  zeros(nnz(gridall),1);
    for nn = 1:15
        opts = detectImportOptions(['inds_',param_list{nn},'.csv']);
        Xv = readtable(['inds_',param_list{nn},'.csv'],opts);
        X((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.tStat;
        pV((nn-1)*1000+1:(nn-1)*1000+size(Xv,1)) = Xv.pValue;
    end
    Te{freq}(gridall==1) = X;
    
    fprintf(sprintf('Loaded frequency %.0f-%.0fHz\n',freql(freq,1),freql(freq,2)))
    
end

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
% close all
freq = 5;
% Te{freq}(Te{freq}==0) = 1;
% sourceant.pow = log10(Te{freq});
for p = 3%1:3
    sourceant =[];

    switch p
        case 1
            sourceant.pow = Tmood{freq};
            fit_parameter = 'mood';
        case 2
            sourceant.pow = Te{freq};
            fit_parameter = 'E';
        case 3
            sourceant.pow = Trpe{freq};
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


crang = [-4 4];
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
    fit_parameter,freql(freq,1),freql(freq,2),max(abs(sourceant.pow(:)))))

% saveas(gcf,sprintf('~/matlab/figures/ELTAgamma_peak5.png'))
end

%%  Plot Linear mixed effects model for sensors
param_list{1} = '001';

if ~exist('sensall','var')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/sensors.mat')
end

data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/';
ff = 5;
 datapath = sprintf('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/pre_mood/powersens_%.0f-%.0fHz/',...
        freql(ff,1),freql(ff,2));
   

    freq = sprintf('powersens_%.0f-%.0fHz',...
        freql(ff,1),freql(ff,2));
    meg_data_name = sprintf('%s.txt',freq);
    meg = dlmread([data_path,meg_data_name]);

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
T.label = sensall;
T.time{1} = 300;
T.sampleinfo = [1 1];
figure; clf; set(gcf,'color','w','position',[176 348 1262 385])

for ii = 1:3%length(fit_parameters)
subplot(1,3,ii)
T.trial{1} = Tfit{ii}; T.avg = Tfit{ii};
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-4 4];
ft_topoplotER(cfg, T)

titlename = fit_parameters{ii};
k = strfind(titlename,'_');
titlename(k) = ' ';
title(titlename)

end

saveas(gcf,sprintf('~/matlab/figures/%s.png',freq))

figure(ff+2); clf
T.trial{1} = mean(megz,2); T.avg =T.trial{1};
cfg = [];
cfg.channel = sensall;
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-.4 .4];
ft_topoplotER(cfg, T)
saveas(gcf,sprintf('~/matlab/figures/%s_avg.png',freq))

% Plot p-values with multiple comparison correction

figure; set(gcf,'color','w')

ii = 2

p = sort(pfit{ii});
semilogy(p)
hold on
semilogy(0.05./(length(sensall):-1:1))
%     grid on
xlabel('sensors')
ylabel('p-value')
legend('p-values','FDR','location','best')
N = nnz(p'<0.05./(length(sensall):-1:1));

title(sprintf('%s: p-value of %.0f sensors < 0.05 (FDR)',fit_parameters{ii},N))

saveas(gcf,sprintf('~/matlab/figures/%s_pvalue.png',freq))


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