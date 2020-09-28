clear all
% close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];
gridres = 5;

%% Load MNI brain and source array
clearvars -except subn gridres

ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));

spmpath = '/spin1/home/linux/liuzzil2/fieldtrip-20190812/external/spm8';

% mri = ft_read_mri([spmpath,'/templates/T1.nii']);
mri = ft_read_mri('~/MNI152_T1_2009c.nii');
nlocs = length(sourcemodel.inside);

%% Combine time courses from all subjects

% Number of datasets per subject
nsets = [1,1,1,1,1,1,2,1,1,1,2,1,1,1,2,2];
% subjects analyzed
slist = [1:9,11,12,14:16]; 

slist = [1,2,4:9,11,14:16]; 


VEpos = cell(sum(nsets(slist)),nlocs);
VEneg = cell(sum(nsets(slist)),nlocs);
ii = 0;
npos = zeros(sum(nsets(slist)),1);
nneg = zeros(sum(nsets(slist)),1);

for sn = slist  %already ran n.6   %Subjects showing enough variation in mood
    
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    cd(data_path)
    name_list = dir;
    data_names = cell(1);
    for k = 1:length(name_list)
        if strncmp(name_list(k).name, data_name, 18) && ~strcmp(name_list(k).name, '24201MMI_mmi3_proc1.ds')
            ii = ii+1;
%             load([data_path,name_list(k).name,'/beamforming/VE_erf_mu1.mat'])
            load([data_path,name_list(k).name,'/beamforming/VE_erf.mat'])
            VEpos(ii,:) = VEp;
            VEneg(ii,:) = VEn;
            npos(ii) = nnz(Rwin);
            nneg(ii) = nnz(~Rwin);
            clear VEp VEn Rwin
        end
    end

end


%% Concatenate data

VEall = cell(1,nlocs);
posarray = false(1,nlocs); 
for ii = 1:length(sourcemodel.inside)
    a = false(1,size(VEpos,1));
    for jj = 1:size(VEpos,1)
        a(jj) = size(VEpos{jj,ii},1) > 1;
    end
    if all(a)
        % zscore in -.2 to 1s window per subject
%         b = cell2mat(VEpos(:,ii)').*npos' +  cell2mat(VEneg(:,ii)').*nneg';         
%         b = abs(zscore(b)); % added this on 07/01/2020
%         bz = zeros(size(b,1),length(slist));
%         zz = 1;
%         for z = 1:length(slist)           
%             if nsets(slist(z))>1
%                 bz(:,z) = ( b(:,zz)*(npos(zz)+nneg(zz)) + b(:,zz+1)*(npos(zz+1)+nneg(zz+1)))/ (sum(npos(zz+[0,1]))+sum(nneg(zz+[0,1]))) ;
%             else
%                 bz(:,z) = b(:,zz);
%             end
%             zz = zz+nsets(z);
%         end
%         VEall{ii} = mean(bz,2);
        
        % Average activity over subjects
        b = cell2mat(VEpos(:,ii)').*npos' +  cell2mat(VEneg(:,ii)').*nneg';         
        b = abs(b);
        bz = zeros(size(b,1),length(slist));
        zz = 1;
        for z = 1:length(slist)           
            if nsets(slist(z))>1
                bz(:,z) = ( b(:,zz) + b(:,zz+1))/ (sum(npos(zz+[0,1]))+sum(nneg(zz+[0,1]))) ;
            else
                bz(:,z) = b(:,zz)/ (npos(zz)+nneg(zz));
            end
            zz = zz+nsets(z);
        end
        VEall{ii} = mean(bz,2);
%         
        % Average activity over all trials
%         b = cell2mat(VEpos(:,ii)').*npos' +  cell2mat(VEneg(:,ii)').*nneg'; 
%         VEall{ii} = sum(abs(b),2)/(sum(npos)+sum(nneg));        
        
        posarray(ii)  = true;    
    end
end


%%

time = linspace(-0.2, 1, 1440 ); % For m = 4*min
% time = linspace(-1,2,3*1200); % For m = 1%max

tstep = 0.05;

for twind = 0.175 %[0.125,0.175,0.225, 0.375]% [0.125,0.175,0.225,0.325, 0.375] %0.075:tstep:0.425  %[0.125,0.175,0.225,0.375]  %0.075:tstep:0.425 

VEerf = zeros(size(VEall));
tt = time>twind & time<(twind+tstep) ; 
for ii = 1:length(posarray)
    if posarray(ii)
        VEerf(ii) = mean(VEall{ii}(tt));
%         VE = VEall{ii};
%         xx = (VE-mean(VE)) > (3*std(VE)); 
%         p = mean(VE(tt & xx'));
%         if isnan(p)
%             VEerf(ii) = 0;
%         else
%             VEerf(ii) = p;
%         end
        
    end
end

% Plot Beamformer
sourceant =[];
sourceant.pow = VEerf;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);

crang = [1 1.8];
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourceout_Int);

set(gcf,'name',sprintf('Win ERF, %.0f-%.0fms',twind*1000,(twind+tstep)*1000))
clc
end


%%
% VE = zeros(size(VEall));
% 
% for ii = 1:length(posarray)
%     if posarray(ii)
%         VEerf = cell2mat(VEpos(:,ii)').*npos' +  cell2mat(VEneg(:,ii)').*nneg';     
%         
%         VEerf = zscore(VEerf);
%         VE(ii) = any(mean(abs(VEerf),2)>1.5);
% %         VE(ii) = max(mean( abs(VEerf),2));
%     end
% end
% 
% sourceant =[];
% sourceant.pow = VE;
% sourceant.dim = sourcemodel.dim;
% sourceant.inside = sourcemodel.inside;
% sourceant.pos = sourcemodel.pos;
% cfg = [];
% cfg.parameter = 'pow';
% sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);
% 
% crang = [0.5 1]; % [1 2]
% cfg = [];
% cfg.method        = 'ortho';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourceout_Int);

%%

x = [27 -75 15
    32 -10 10
    37 10 15
    2 -19 47 
    2 48 1];

labelst = {'Left occipital';'left posterior insula';'left frontal insula';'MCC';'ACC arb.'};

% 
% x = [
%     -4 23 30
%     -15 -5 6    
%     30 -46 46
%     38 -56 0
%     ];
% 
% labelst = {'ACC';'Right striatum';'lS';'posterior temporal'};



% 
% x = [27 -75 15
%     32 -10 10
%     37 10 15
%     -4 23 30
%     -15 -5 6
%     2 -19 47
%     30 -46 46
%     38 -56 0
%     2 48 1];
% 
% labelst = {'Left occipital';'left posterior insula';'left frontal insula';'ACC';'Right striatum';'MCC';'lS';'posterior temporal';'ACC arb.'};





newf = 250;
% timed = linspace(-1,2,newf*3);
timed = linspace(-0.2,1,newf*1.2);
voxpos = round(sourcemodel.pos*10); % convert to mm
figure; clf
for ii = 1:size(x,1)
    
    x1 = (abs(voxpos(:,1) - x(ii,1))<=5);
    x2 = (abs(voxpos(:,2) - x(ii,2))<=5);
    x3 = (abs(voxpos(:,3) - x(ii,3))<=5);
    xx = x1 & x2 & x3;
    
    %
    %
    %     b = cell2mat(VEpos(:,xx)').*npos' +  cell2mat(VEneg(:,xx)').*nneg';
    %     b = reshape(b, length(time), length(find(xx)));
    %
    %     s = corr(b);
    %     b = b.*sign(s(1,:));  
    %     VE = mean(b,2);
    %     VE = zscore(VE);
    b = cell2mat(VEpos(:,xx)').*npos' +  cell2mat(VEneg(:,xx)').*nneg';
    b = reshape(b, length(time), length(find(xx)),size(b,2));
    b = zscore(b);
    
    bzm = zeros(size(b,1),size(b,3));
    for z = 1:size(b,3)
        bz = b(:,:,z);
        s = corr(bz);
        bz = bz.*sign(s(1,:));
        bzm(:,z) = squeeze(mean(bz,2));
    end
    
    
    bz = zeros(size(b,1),length(slist));
    zz = 1;
    for z = 1:length(slist)
        if nsets(slist(z))>1
            bz(:,z) = ( bzm(:,zz)*(npos(zz)+nneg(zz)) + bzm(:,zz+1)*(npos(zz+1)+nneg(zz+1)))/ (sum(npos(zz+[0,1]))+sum(nneg(zz+[0,1]))) ;
        else
            bz(:,z) = bzm(:,zz);
        end
        zz = zz+nsets(z);
    end
    
    bz = resample(bz,newf,1200);
    s = corr(bz);
    bz = bz.*sign(s(1,:));
    
    VE = mean(bz,2);
    
       
    if mean(VE(timed>0.2 & timed < 0.35))> 0%abs(max(VE)) > abs(min(VE))
        plot(timed*1e3,VE)
    else
        plot(timed*1e3,-VE)
    end
    
    hold on
end

xlabel('Time (ms)'); ylabel('zscore'); %ylabel('nAm')
legend(labelst)
grid on
% xlim([0 500])


%% Positive and Negative RPE separately
% x = [27 -75 15
%     32 -10 10
%     37 10 15
%     -4 23 30
%     -15 -5 6
%     2 -19 47
%     30 -46 46
%     38 -56 0
%     2 48 1];
% 
% labelst = {'Left occipital';'left posterior insula';'left frontal insula';'ACC';'Right striatum';'MCC';'lS';'posterior temporal';'ACC arb.'};

x = [28 5 15
    10 -47 5
    9 -26 42
    41 -13 10
    27 -80 15
    -24 1 20
    9 -50 -21
    4 18 -11
    2 32 17
    ];
labelst = {'Left striatum';'PCC';'medial cingulate';'Left insula';'Left Occipital';'Right striatum';'Cerebellum';'Nucleus Acumbens';'ACC'};

x = [27 -75 15
    32 -10 10
    37 10 15
    2 -19 47 
    2 48 1];

labelst = {'Left occipital';'left posterior insula';'left frontal insula';'MCC';'ACC arb.'};


voxpos = round(sourcemodel.pos*10); % convert to mm

for ii = 1:size(x,1)
    
    x1 = (abs(voxpos(:,1) - x(ii,1))<=5);
    x2 = (abs(voxpos(:,2) - x(ii,2))<=5);
    x3 = (abs(voxpos(:,3) - x(ii,3))<=5);
    xx = x1 & x2 & x3;
    
    % Positive RPE
    b = cell2mat(VEpos(:,xx)');
    b = reshape(b, length(time), length(find(xx)),size(b,2));
    b = zscore(b);
    
    bzm = zeros(size(b,1),size(b,3));
    for z = 1:size(b,3)
        bz = b(:,:,z);
        s = corr(bz);
        bz = bz.*sign(s(1,:));
        bzm(:,z) = squeeze(mean(bz,2));
    end
    
    
    bz = zeros(size(b,1),length(slist));
    zz = 1;
    for z = 1:length(slist)
        if nsets(slist(z))>1
            bz(:,z) = ( bzm(:,zz)*(npos(zz)) + bzm(:,zz+1)*(npos(zz+1)))/ (sum(npos(zz+[0,1]))) ;
        else
            bz(:,z) = bzm(:,zz);
        end
        zz = zz+nsets(z);
    end
    
    bz1 = resample(bz,newf,1200);
    s = corr(bz1);
    bz1 = bz1.*sign(s(1,:));
    
    VE1 = mean(bz1,2);
%     VE1 = bz(:,sn);
    
    
    % Negative RPE
    b = cell2mat(VEneg(:,xx)');
    b = reshape(b, length(time), length(find(xx)),size(b,2));
    b = zscore(b);
    
    bzm = zeros(size(b,1),size(b,3));
    for z = 1:size(b,3)
        bz = b(:,:,z);
        s = corr(bz);
        bz = bz.*sign(s(1,:));
        bzm(:,z) = squeeze(mean(bz,2));
    end
    
    
    bz = zeros(size(b,1),length(slist));
    zz = 1;
    for z = 1:length(slist)
        if nsets(slist(z))>1
            bz(:,z) = ( bzm(:,zz)*(nneg(zz)) + bzm(:,zz+1)*(nneg(zz+1)))/ (sum(nneg(zz+[0,1]))) ;
        else
            bz(:,z) = bzm(:,zz);
        end
        zz = zz+nsets(z);
    end
    bz2 = resample(bz,newf,1200);
    s = corr(bz2);
    bz2 = bz2.*sign(s(1,:));
    
    VE2 = mean(bz2,2);
%     VE2 = bz(:,sn);
    
    figure(ii+5); clf
    
    plot(timed*1e3,VE1)
    hold all 
    if corr(VE1,VE2)< 0%abs(max(VE)) > abs(min(VE))        
        VE2 = -VE2;
    end
    plot(timed*1e3,VE2)
    
    fill([timed fliplr(timed)]'*1e3,[VE1-std(bz1,0,2)/sqrt(length(slist));...
        flipud(VE1+std(bz1,0,2)/sqrt(length(slist)))],...
        [0 0 1],'facealpha',0.2,'edgecolor','none')
    
    fill([timed fliplr(timed)]'*1e3,[VE2-std(bz2,0,2)/sqrt(length(slist));...
        flipud(VE2+std(bz2,0,2)/sqrt(length(slist)))],...
         [1 0 0],'facealpha',0.2,'edgecolor','none')
    
    title(labelst{ii})
    xlabel('Time (ms)'); ylabel('zscore'); %ylabel('nAm')
    legend('Positive RPE','Negative RPE')
    grid on
    set(gcf,'color','white')
end






%% PCA?
VE = zeros(length(time),nnz(posarray));

jj  =0;
for ii = 1:length(posarray)
    if posarray(ii)
    jj = jj+1;
    b = (cell2mat(VEpos(:,ii)').*npos' +  cell2mat(VEneg(:,ii)').*nneg') ./ (npos+nneg)';
%     b = cell2mat(VEneg(:,ii)');
%     b = cell2mat(VEpos(:,ii)') -  cell2mat(VEneg(:,ii)');
    
%     b = zscore(b);
   
    s = corr(b);
    b = b.*sign(s(1,:));
        
    
    bz = zeros(size(b,1),length(slist));
    zz = 1;
    for z = 1:length(slist)
        if nsets(slist(z))>1
            bz(:,z) = ( b(:,zz)*(npos(zz)+nneg(zz)) + b(:,zz+1)*(npos(zz+1)+nneg(zz+1)))/ (sum(npos(zz+[0,1]))+sum(nneg(zz+[0,1]))) ;
%             bz(:,z) = ( b(:,zz)*(npos(zz)) + b(:,zz+1)*(npos(zz+1)))/ (sum(npos(zz+[0,1]))) ;
%             bz(:,z) = ( b(:,zz)*(nneg(zz)) + b(:,zz+1)*(nneg(zz+1)))/ (sum(nneg(zz+[0,1]))) ;

        else
            bz(:,z) = b(:,zz);
          
        end
        zz = zz+nsets(z);
    end
    
%     s = corr(bz);
%     bz = bz.*sign(s(1,:));
    
    VE(:,jj) = mean(bz,2);
    end
 
end

[coeff,score] = pca(VE);

%%
figure; plot(time,score(:,1:3))
title('PCA components')


for ii =1:3
coeffplot = zeros(size(posarray));
coeffplot(posarray) = coeff(:,ii);

sourceant =[];
sourceant.pow = abs(coeffplot);
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri);

crang = []; 
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;

ft_sourceplot(cfg, sourceout_Int);
title(sprintf('Principal compoent %d',ii))
end