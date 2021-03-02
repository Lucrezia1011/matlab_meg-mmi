% Based on mmi_LTA_aal_prep_RUN
clear all
close all
clc
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
subexclude = [10,24];

roiopt = 'AAL'; 
switch roiopt
    case 'AAL'
        subexclude = [subexclude,26,49,53];
end

Nlist(subexclude) = []; 
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

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

twind = [-.2,1];

for ii = 1:length(data_list) 
    data_name = data_list{ii};
    sub = data_name(5:9);
    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
%     mmiAALprep(data_list{ii},twind,roiopt)

    gridres=8; % voxel distance in mm
    mu = 0.01; % covariance regularization parameter: fraction of maximum singular value.
    sampfreq = 50; % downsampled frequency, 50Hz
    
    save_name = sprintf('%s/evoked_outcome_mu-%.2f_res-%.0fmm_fs-%.0fHz',...
        processing_folder,mu,gridres,sampfreq);
    if exist(save_name,'file')
        eval(sprintf('!mv %s %s',save_name,[save_name,'.mat'])) 
    elseif ~exist([save_name,'.mat'],'file')
        mmiEvokedPrep(data_name,twind,gridres,mu,sampfreq)
    end
end
return
%%
freq = 'outcome';% 'outcome', 'cue' , 'choice'
% save(save_name ,'Ybeta','Ytheta','Y','ltvall')

if strcmp(roiopt,'sens')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    cd('/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/')
    hdr = ft_read_header(data_list{1});
    channelsall = hdr.label(strcmp(hdr.chantype,'meggrad'));
end
time = linspace(twind(1),twind(2),60);

if ~exist('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/mni_grid.txt','file')
    for sn = 1:length(data_list) % all subjects with continuos recordings and latent variables

        data_name = data_list{sn};   
        sub = data_name(5:9);
        processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        load([processing_folder,'leadfields_',num2str(gridres),'mm.mat'])
        if sn ==1
            gridall = grid.inside;
        else
            gridall = gridall & grid.inside;

        end
    end
    dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/mni_grid.txt',gridall)
else
    gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/mni_grid.txt');
end

%%
r = 0;
Yall = [];
ltv = [];

Ym = [];

ntrials = [];
ntrialsOut = [];

for sn = 1:length(data_list) % all subjects with continuos recordings and latent variables
    
    data_name = data_list{sn};   
    sub = data_name(5:9);
    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save_name = [processing_folder,'evoked_',freq,'_',roiopt,'.mat'];
    load(save_name)
    
    save_name = sprintf('%s/evoked_outcome_mu-%.2f_res-%.0fmm_fs-%.0fHz.mat',...
        processing_folder,mu,gridres,sampfreq);
    load(save_name)
    
    load([processing_folder,'leadfields_',num2str(gridres),'mm.mat'])
    Yout = zeros(length(gridall),size(VEall,2),size(VEall,3));
    Yout(grid.inside==1,:,:) = VEall;
    Yout = Yout(gridall==1,:,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = r+1;
    ltvout.recording = repmat(r,size(ltvout,1),1);
   
     if strcmp(roiopt,'sens')
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

        [~,~,ind]= intersect(channels,channelSub);
        Yout = Yout(ind,:,:);
        
        % Can only use this for average response (over subject or time)
        hdr = ft_read_header(data_name);
        cfg_neighb        = [];
        cfg_neighb.method = 'template';%'distance';
        cfg_neighb.channel = channelSub;
        cfg_neighb.template = 'CTF275_neighb.mat';
        neighbours        = ft_prepare_neighbours(cfg_neighb, hdr);
        
        ntrials = size(Yout,3);
        
        
     end
    
    
    Yall = cat(3,Yall,Yout);
    % Check sign of evoked response
%     Ym = cat(3,Ym, cat(2,mean(Ycue,3),mean(Yout,3),mean(Ychoice,3)));
    Ym = cat(3,Ym, cat(2,mean(Yout,3)));
    ntrialsOut = cat(1,ntrialsOut,size(Yout,3));
    
    if isempty(ltv)
        ltv = ltvout;
    else
        ltv(end+(1:size(ltvout,1)),:) = ltvout;
    end
    
    fprintf('Done %.0f/%.0f\n',sn,length(data_list))
    
    
end


Ym0 = Ym;
% ltv = flipud(ltv);
clear ltvchoice ltvcue ltvout Y Yalpha Ybeta Ytheta Ycue Yout Ychoice
%% Sign flip for source localized evoked responses

% Yu = squeeze(Ym(ii,:,:))';
Yu = squeeze(sum(Ym,1))';

Yu = zscore(Yu,0,2);

%    % first pass
Yav = [];
C = corr(Yu');
C = triu(C,1);
C( abs(C-1) < 0.01) = 0;
[~,ind] = sort(abs(C(:)),'descend');
m = sign(C(ind(1)));
[i,j] = ind2sub(size(C),ind(1));
Yav = cat(1,Yav, Yu(i,:), m*Yu(j,:));
s = zeros(size(C,1),1);
s(i) = 1;
s(j) = m;

while nnz(s) < size(Yu,1)
    C = corr(mean(Yav,1)',Yu');
    [~,ind] = sort(abs(C(:)),'descend');
    inds = find(s);
    z = 1;
    while any(ind(z) == inds)
        z = z+1;
    end
    
    m = sign(C(ind(z)));
    s(ind(z)) = m;
    Yav =  cat(1,Yav, m*Yu(ind(z),:));
end

Y0 = Yall;
zzo = 0;

Sout = zeros(size(Yall));

for z = 1:size(Yu,1)
    Sout(:,:,zzo+(1:ntrialsOut(z))) = s(z);
    zzo = zzo + ntrialsOut(z);
end
Y0 = Y0.*Sout;

figure; subplot(121);
imagesc(mean(Ym,3)); caxis([-1 1]*1e-13)
subplot(122);
imagesc(mean(Y0,3)); caxis([-1 1]*1e-13)

%% Write data
N= 10; % save as N files
nstep = ceil(size(Y0,1)/N);
for n = 1:N
    npoints = ((n-1)*nstep +1) : n*nstep; 
    npoints(npoints>size(Y0,1)) =[];
    Yall = permute(Y0(npoints,:,:),[2,1,3]);
    Yall = reshape(Yall,size(Yall,1)*size(Yall,2),size(Yall,3)); % nrois * npoints * ntrials

    dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/meg_trials_evoked_',freq,'_',num2str(n),'.txt'],Yall);
end
writetable(ltv,['/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/latent_vars_evoked_',freq,'.csv']);


%% Read data
Yall = [];
N= 10; % save as N files
for n = 1:N

    Y = dlmread(['/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/lme_RPE_sum/voxs_',num2str(n),'.txt']);
    Y = reshape(Y,[60,length(Y)/60]);
    Yall = cat(2,Yall,Y);
end
gridall = dlmread('/data/MBDU/MEG_MMI3/results/mmiTrial_grid/evoked_outcome/mni_grid.txt');

mri_mni = ft_read_mri('~/fieldtrip-20190812/external/spm8/templates/T1.nii','dataformat','nifti');
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';

load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
sourcemodel.coordsys = 'mni';

clusternull2 = zeros(size(gridall,1),60);
clusternull2(gridall==1,:) = Yall';
img = reshape(clusternull2,[sourcemodel.dim,60]);


E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1;

[tfced] = matlab_tfce_transform_gridtime(img,H,E,dh);
tfced = tfced - matlab_tfce_transform_gridtime(-img,H,E,dh);

%%
close all

% [s,ind] = sort(abs(tfced(:)),'descend');
% [~,~,~,tt] = ind2sub(size(img),ind(1:10));
[m,indm]=max(tfced(:));
[~,~,~,ttm] = ind2sub(size(img),indm);

[m,indn]=min(tfced(:));
[~,~,~,ttn] = ind2sub(size(img),indn);

tt=  [ttm,ttn]
time = linspace(-200,1000,60);
% for t =[650]
% tt = t+[0 50] ;
% [~,ii]= min(abs(time - tt'),[],2);
% 
% sourceant.pow  = mean(tfced(:,:,:,ii(1):ii(2)),4);
% T  = img;
% T(abs(tfced)<6) =0;
% sourceant.pow  = mean(T(:,:,:,ii(1):ii(2)),4);

for t = 1:2

% sourceant.pow  =tfced(:,:,:,tt(t));
T  = img;
T(abs(tfced)<3) =0;
sourceant.pow  = T(:,:,:,tt(t));

sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;


cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceant , mri_mni);
sourceout_Int.pow(~sourceout_Int.inside) = 0;
sourceout_Int.coordsys = 'mni';


crang = [-5 5];
% crang = [thresh max(sourceant.pow)];
cfg = [];
cfg.method        = 'ortho'; %'ortho'
if max(sourceout_Int.pow(:)) > -min(sourceout_Int.pow(:))
    cfg.location   = 'max';
else
    cfg.location   = 'min';
end
%     cfg.location = [-2 -40 30]; % peak for mu =0.05
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
cfg.atlas = '~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii';

ft_sourceplot(cfg, sourceout_Int);
% title(sprintf('Time %.0f-%.0fms',time(ii(1)),time(ii(2))))
title(sprintf('Time %.0fms',time(tt(t))))

cfg.method        = 'slice'; %'ortho'
ft_sourceplot(cfg, sourceout_Int);
% title(sprintf('Time %.0f-%.0fms',time(ii(1)),time(ii(2))))
title(sprintf('Time %.0fms',time(tt(t))))

end