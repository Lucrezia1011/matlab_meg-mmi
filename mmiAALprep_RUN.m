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
Nlist([10,24,26,49,53]) = [];
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
% 
% roiopt = 'g';
% for ii = 30:length(data_list) 
%     data_name = data_list{ii};
%     sub = data_name(5:9);
%     processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
%     save_name = [processing_folder,'evoked_outcome_AAL.mat'];
%     mmiAALprep(data_list{ii},twind)
% end
% return
%%
freq = 'outcome';
% save(save_name ,'Ybeta','Ytheta','Y','ltvall')
r = 0;
Yall = [];
YOut = [];

ltv = [];
ltvOut = [];

Ym = [];

ntrials = [];
ntrialsOut = [];
for sn = 1:length(data_list) % all subjects with continuos recordings and latent variables
    
    data_name = data_list{sn};   
    sub = data_name(5:9);
    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save_name = [processing_folder,'evoked_',freq,'_AAL.mat'];
    load(save_name)
    
    switch freq
        case 'cue'
            ltvout = ltvcue;
            Yout = Ycue;
        case 'choice'
            ltvout = ltvchoice;
            Yout = Ychoice;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = r+1;
    ltvout.recording = repmat(r,size(ltvout,1),1);
   
    YOut = cat(3,YOut,Yout);
    % Check sign of evoked response
%     Ym = cat(3,Ym, cat(2,mean(Ycue,3),mean(Yout,3),mean(Ychoice,3)));
    Ym = cat(3,Ym, cat(2,mean(Yout,3)));
    ntrialsOut = cat(1,ntrialsOut,size(Yout,3));
    
    if isempty(ltvOut)
        ltvOut = ltvout;
    else
        ltvOut(end+(1:size(ltvout,1)),:) = ltvout;
    end
    
    
    
end

Ym0 = Ym;
% ltv = flipud(ltv);
clear ltvchoice ltvcue ltvout Y Yalpha Ybeta Ytheta Ycue Yout Ychoice
%% Sign flip for source localized evoked responses
% S = zeros(size(Yall));

Sout = zeros(size(YOut));
% Scue = zeros(size(YCue));
% Schoice = zeros(size(YChoice));

% Ym = Ym0(:,[100:300,390+(80:200)],:);

for ii = 1:size(Ym,1)
    Yu = squeeze(Ym(ii,:,:))';
%     [Yu,IA,IC] = unique(Y','rows');
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
    
    
    zzo = 0;
    zzcu = 0;
    zzch = 0;
    
    for z = 1:size(Yu,1)
        Sout(ii,:,zzo+(1:ntrialsOut(z))) = s(z);
        zzo = zzo + ntrialsOut(z);
%         Scue(ii,:,zzcu+(1:ntrialsCue(z))) = s(z);
%         zzcu = zzcu + ntrialsCue(z);
%         Schoice(ii,:,zzch+(1:ntrialsChoice(z))) = s(z);
%         zzch = zzch + ntrialsChoice(z);
    end
    
    
end

% Yall = Yall.*S;


figure; set(gcf,'color','w')
% subplot(231); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
% subplot(232); imagesc(mean(YChoice,3)); caxis([-1 1]*4e-12)
subplot(131); imagesc(mean(YOut,3)); caxis([-1 1]*.3)

% YCue = YCue.*Scue;
% YChoice = YChoice.*Schoice;
YOut = YOut.*Sout;

% subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
% subplot(235); imagesc(mean(YChoice,3)); caxis([-1 1]*4e-12)
subplot(132); imagesc(mean(YOut,3)); caxis([-1 1]*.3)


% Flip ROIs to match each other
% Yu = cat(2,mean(YCue,3),mean(YChoice,3),mean(YOut,3));
Yu = mean(YOut,3);

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
npoints = size(Sout,2);
Sout = repmat(s,[1,npoints,sum(ntrialsOut)]);
% npoints = size(Scue,2);
% Scue = repmat(s,[1,npoints,sum(ntrialsCue)]);
% npoints = size(Schoice,2);
% Schoice = repmat(s,[1,npoints,sum(ntrialsChoice)]);

% YCue = YCue.*Scue;
% YChoice = YChoice.*Schoice;
YOut = YOut.*Sout;

% subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
% subplot(235); imagesc(mean(YChoice,3));  caxis([-1 1]*4e-12)
subplot(133); imagesc(mean(YOut,3));  caxis([-1 1]*.3)



%% PCA evoked responses experiment
% if ~exist('atlas','var')    
%     atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
% end
% 
% 
% n = 55;
% time =  linspace(-3,3,600);
% Y = squeeze(YAlpha(n,:,:));
% [coeff,score] = pca(Y);
% 
% figure(1); clf 
% subplot(211)
% plot(time,mean(Y,2))
% title(atlas.tissuelabel{n})
% subplot(212)
% plot(time,score(:,1))
% hold all
% plot(time,score(:,2))
% % plot(time,score(:,3))
% % plot(time,score(:,4))
% title('PCA')
%%

YOut = permute(YOut,[2,1,3]);
YOut = reshape(YOut,size(YOut,1)*size(YOut,2),size(YOut,3)); % nrois * npoints * ntrials

% YCue = permute(YCue,[2,1,3]);
% YCue = reshape(YCue,size(YCue,1)*size(YCue,2),size(YCue,3)); % nrois * npoints * ntrials
% 
% YChoice = permute(YChoice,[2,1,3]);
% YChoice = reshape(YChoice,size(YChoice,1)*size(YChoice,2),size(YChoice,3)); % nrois * npoints * ntrials


%% Write data

dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/meg_trials_evoked_',freq,'.txt'],YOut);
% dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_',freq,'_choice.txt'],YChoice);
% dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_',freq,'_cue.txt'],YCue);

writetable(ltvOut,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/latent_vars_evoked_',freq,'.csv']);
% writetable(ltvChoice,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_',freq,'_choice.csv']);
% writetable(ltvCue,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_',freq,'_cue.csv']);




