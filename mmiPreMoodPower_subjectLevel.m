clear all
close all

% start fieldtrip
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
%% Do with average power
out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood';

ltvMood = readtable([out_path,'/latent_vars.csv']);

freqband=[ 25 40];
mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
% mri_mni = ft_read_mri('~/MNI152_T1_2009c.nii'); % in mni coordinates

ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
gridres = 5;
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/mni_grid.txt');
Pall = dlmread(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
    freqband(1),freqband(2)));
% Pall = dlmread(sprintf('%s/powergrid_%.0f-%.0fHz_mu0.2.txt',out_path,...
%     freqband(1),freqband(2)));

bestfit_name = '/data/MBDU/MEG_MMI3/data/behavioral/closed_LTA_coefs.csv';
opts = detectImportOptions(bestfit_name);
bf_pars = readtable(bestfit_name,opts); 

meginfoname = '~/MEG_participantsinfo.csv';
opts = detectImportOptions(meginfoname);
meginfo = readtable(meginfoname,opts); 

subs = unique(ltvMood.subject);
Ps = zeros(size(Pall,1),length(subs));
rn = cell(1,length(subs));
rnm = zeros(1,length(subs));
w_LTE = zeros(1,length(subs));
w_LTR = zeros(1,length(subs));
M0 = zeros(1,length(subs));
MDD = zeros(1,length(subs));
for s = 1:length(subs)
    Ps(:,s) = mean(Pall(:,ltvMood.subject == subs(s)),2);
    recs = ltvMood.recording(ltvMood.subject == subs(s));
    nrecs = unique(recs);
    for n = 1:length(nrecs)
        rn{s}(n) = nnz(recs == nrecs(n));
    end
    rnm(s) = mean(rn{s});
    bestfit_sub = bf_pars(bf_pars.Var1 == subs(s),:);
    w_LTE(s)  = bestfit_sub.w_LTE;
    w_LTR(s)  = bestfit_sub.w_LTR;
    M0(s)  = bestfit_sub.m_0;
    MDD(s) = strcmp(meginfo.Group(meginfo.SDAN==subs(s)),'MDD');
end

r = zeros(1,size(Ps,1));
pv= zeros(1,size(Ps,1));
for n = 1:size(Ps,1)
y =  Ps(n,:);
x  = w_LTE;
% p = polyfit(x,y,1);
% yp = polyval(p,x);
% yres = y-yp;
% SSresid = sum(yres.^2);
% %Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
% SStotal = (length(yp)-1) * var(y);
% % Compute R2 using the formula given in the introduction of this topic:
% rsq = 1 - SSresid/SStotal;
[r(n),pv(n)]=corr(x',y');
end

figure;clf;set(gcf,'color','w')

subplot(221)
scatter(M0,w_LTE)
p = polyfit(M0,w_LTE,1);
hold on
plot(M0,M0*p(1)+p(2))
xlabel('M_0'); ylabel('\beta_E')
title(sprintf('m=%.3f, r=%.1f',p(1),corr(M0',w_LTE')))

subplot(222)
scatter(w_LTR,w_LTE)
p = polyfit(w_LTR,w_LTE,1);
hold on
plot(w_LTR,w_LTR*p(1)+p(2))
xlabel('\beta_R'); ylabel('\beta_E')
title(sprintf('m=%.2f, r=%.1f',p(1),corr(w_LTR',w_LTE')))

subplot(223)
boxplot(M0,MDD)
set(gca,'XTickLabel',{'HV';'MDD'})
ylabel('M_0')
title(sprintf('M_0 of HVs and MDDs are\nnot signtifcantly different'));

subplot(224)
boxplot(w_LTE,MDD)
set(gca,'XTickLabel',{'HV';'MDD'})
ylabel('\beta_E')
title(sprintf('\\beta_E of HVs and MDDs are\nnot signtifcantly different'));

%% 
% cortical mask
aal = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d5mm']));
aal = ft_convert_units(aal,sourcemodel.unit);
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodelAAL = ft_sourceinterpolate(cfg, aal, sourcemodel);

% Plot Correlation

R = zeros(size(gridall));
R(gridall==1) =r;
sourceant.pow  = R;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';

crang = [0.45 0.65];
% crang = [min(sourceant.pow(sourceant.pow>0)) max(sourceant.pow(:))];
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

% X = reshape(sourcemodel.pos(:,1),sourcemodel.dim);
% Y = reshape(sourcemodel.pos(:,2),sourcemodel.dim);
% Z = reshape(sourcemodel.pos(:,3),sourcemodel.dim);
% fv = patch(isosurface(sourcemodelAAL.tissue>0 & sourcemodelAAL.tissue<=90,0.5));
% isonormals(sourcemodelAAL.tissue,fv)
% fv.FaceColor = [0.5 0.5 0.5];
% fv.EdgeColor = 'none';
% daspect([1 1 1])
% view(3); 
% axis tight
% camlight 
% lighting gouraud

[X,Y,Z] = meshgrid(1:aal.dim(2),1:aal.dim(1),1:aal.dim(3));
fv = (isosurface(X,Y,Z,aal.tissue>0 & aal.tissue<=90,0.5));
IN = inpolyhedron(fv,[X(:),Y(:),Z(:)]); %QPTS = Nx3
% Plot the mask
sourceout_Int = mri_mni;
sourceout_Int.pow  = aal.tissue;
% sourceout_Int.pow(aal.tissue>90 | aal.tissue<1) = 0;
sourceout_Int.pow(aal.tissue==0) = 100;
sourceout_Int.pow(~IN) = 0;
sourceout_Int.coordsys = 'mni';
sourceout_Int.inside = aal.tissue<=90 & aal.tissue>=1;
crang = [];
% crang = [min(sourceant.pow(sourceant.pow>0)) max(sourceant.pow(:))];
cfg = [];
cfg.method        = 'ortho'; %'ortho'
cfg.funparameter = 'pow';
cfg.maskparameter = 'inside';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);

%% TFCE with random permutations, AAL only
E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1; 
C =26;
R = zeros(size(gridall));
R(gridall==1) =r;
R(sourcemodelAAL.tissue>90  | sourcemodelAAL.tissue==0) = 0;
R = reshape(R,sourcemodel.dim);
tfce= matlab_tfce_transform(R,H,E,C,dh); % C=26 default setting
tfcen= matlab_tfce_transform(-R,H,E,C,dh); % C=26 default setting
tfce = tfce-tfcen;

Nr = 5000; % number of random iterations
tfce_rand = cell(1,Nr);
parfor nr = 1:Nr
    rn = zeros(1,size(Ps,1));
    x  = w_LTE(randperm(length(w_LTE)))';
    for n = 1:size(Ps,1)
        y =  Ps(n,:);
        rn(n)=corr(x,y');
    end
    Rn = zeros(size(gridall));
    Rn(gridall==1) =rn;
    Rn(sourcemodelAAL.tissue>90 | sourcemodelAAL.tissue==0 ) = 0;
    Rn = reshape(Rn,sourcemodel.dim);
    tfcer= matlab_tfce_transform(Rn,H,E,C,dh);
    tfcen= matlab_tfce_transform(-Rn,H,E,C,dh);
    tfcer = tfcer-tfcen;
    tfce_rand{nr} = max(tfcer(:));
    clc
    fprintf('Done %.0f perc.\n',nr/Nr*100)
end  


tfce_rand = cell2mat(tfce_rand); 
tfce_rand = sort(tfce_rand(:));
% thresh = tfce_rand(0.99*1000*nnz(gridall));
thresh = tfce_rand(0.99*Nr);

%%
sourceant.pow  = R;
sourceant.pow(abs(tfce)<thresh)  = 0; % one tailed
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';

% crang = [];
crang = [min(sourceant.pow(sourceant.pow>0)) max(sourceant.pow(:))];

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




cfg = [];
cfg.parameter = 'pow';
mask = ft_read_mri('/data/MBDU/MEG_MMI3/data/hanna_acc_mask/ACC_MoodQ_LTE_mask_0001_mni.nii');

mask.pow = mask.anatomy; %mask.anatomy>0.9;
cfg = [];
cfg.parameter = 'pow';

sourceout_IntH  = ft_sourceinterpolate(cfg, mask , mri_mni);

sourceout_IntH.pow(~sourceout_IntH.inside) = 0;
sourceout_IntH.coordsys = 'mni';

% ccl = linspace( min(sourceout_Int.pow(sourceout_Int.pow>0)),  max(sourceout_Int.pow(:)),256);
ccl = linspace( 0.4, 0.61,256);


sourceout_Int.pow(sourceout_Int.pow>0 & sourceout_IntH.pow>0) = ccl(end);
%     -sourceout_Int.pow(sourceout_Int.pow>0  & sourceout_IntH.pow>0) ; % make =-1 to plot 
sourceout_Int.pow(sourceout_Int.pow==0 & sourceout_IntH.pow>0) =  ccl(end-1); % make =-1 to plot 

cc = colormap(hot);
cc(end,:) = [0 1 0];
cc(end-1,:) = [0 0 1];

crang = [ccl(1) ccl(end)];
cfg = [];
cfg.method        = 'ortho'; %'slice'
cfg.location   = 'max';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = cc;
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
ft_sourceplot(cfg, sourceout_Int);

saveas(gcf,sprintf('~/matlab/figures/ACCoverlp_sublevel_0.01TFCE_AALonly.tif'))

[~,indP] = min(sum((sourcemodel.pos-[-2 30 0]/10).^2,2));
indG = 1:length(gridall);
indG = indG(gridall==1);
indf = find(indG == indP);
figure; set(gcf,'color','w');
scatter(w_LTE,Ps(indf,:))
pp = polyfit(w_LTE,Ps(indf,:),1);
hold on
plot(w_LTE,w_LTE*pp(1)+pp(2))
xlabel('\beta_E'); ylabel('25-40Hz power (a.u. w''Cw/w''C_{noise}w)')
title(sprintf('Corr. of \\beta_{E} and 25-40Hz power at ACC peak\nr = %.2f, p-value = %.1e',r(indf),pv(indf)))


%% TFCE with random permutations, all voxels
close all
E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1; 
C =26;
R = zeros(size(gridall));
R(gridall==1) =r;
% R(sourcemodelAAL.tissue>90  | sourcemodelAAL.tissue==0) = 0;
R = reshape(R,sourcemodel.dim);
tfce= matlab_tfce_transform(R,H,E,C,dh); % C=26 default setting
tfcen= matlab_tfce_transform(-R,H,E,C,dh); % C=26 default setting
tfce = tfce-tfcen;

Nr = 5000; % number of random iterations
tfce_rand = cell(1,Nr);
parfor nr = 1:Nr
    rn = zeros(1,size(Ps,1));
    x  = w_LTE(randperm(length(w_LTE)))';
    for n = 1:size(Ps,1)
        y =  Ps(n,:);
        rn(n)=corr(x,y');
    end
    Rn = zeros(size(gridall));
    Rn(gridall==1) =rn;
%     Rn(sourcemodelAAL.tissue>90 | sourcemodelAAL.tissue==0 ) = 0;
    Rn = reshape(Rn,sourcemodel.dim);
    tfcer= matlab_tfce_transform(Rn,H,E,C,dh);
%     tfcen= matlab_tfce_transform(-Rn,H,E,C,dh);
%     tfcer = tfcer-tfcen;
    tfce_rand{nr} = max(tfcer(:));
    clc
    fprintf('Done %.0f perc.\n',nr/Nr*100)
end  


tfce_rand = cell2mat(tfce_rand); 
tfce_rand = sort(tfce_rand(:));
% thresh = tfce_rand(0.99*1000*nnz(gridall));
thresh = tfce_rand(0.99*Nr);

%% ACC mask from Hanna's paper

sourceant.pow  = R;
sourceant.pow(tfce<thresh)  = 0; % one tailed
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';


cfg = [];
cfg.parameter = 'pow';
mask = ft_read_mri('/data/MBDU/MEG_MMI3/data/hanna_acc_mask/ACC_MoodQ_LTE_mask_0001_mni.nii');

mask.pow = mask.anatomy; %mask.anatomy>0.9;
cfg = [];
cfg.parameter = 'pow';

sourceout_IntH  = ft_sourceinterpolate(cfg, mask , mri_mni);

sourceout_IntH.pow(~sourceout_IntH.inside) = 0;
sourceout_IntH.coordsys = 'mni';

% ccl = linspace( min(sourceout_Int.pow(sourceout_Int.pow>0)),  max(sourceout_Int.pow(:)),256);
ccl = linspace( 0.2, 0.632,256);


sourceout_Int.pow(sourceout_Int.pow>0 & sourceout_IntH.pow>0) = ccl(end);
%     -sourceout_Int.pow(sourceout_Int.pow>0  & sourceout_IntH.pow>0) ; % make =-1 to plot 
sourceout_Int.pow(sourceout_Int.pow==0 & sourceout_IntH.pow>0) =  ccl(end-1); % make =-1 to plot 

cc = colormap(hot);
cc(end,:) = [0 1 0];
cc(end-1,:) = [0 0 1];

crang = [ccl(1) ccl(end)];
cfg = [];
cfg.method        = 'ortho'; %'slice'
cfg.location   = 'max';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = cc;
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';
ft_sourceplot(cfg, sourceout_Int);

saveas(gcf,sprintf('~/matlab/figures/ACCoverlp_sublevel_0.01TFCE.tif'))