clear all
close all
clc

subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn = 1:16 %[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

% twind = [-.2,1];
% evokedopt = true; % true or false
% inducedopt = {'theta';'alpha';'beta';'delta'}; % cell with freq limitis or empty
% roiopt = 'g';
% for ii = 1:length(param_list) % redo 15
%     mmi_LTA_aal_prep(param_list{ii},twind,evokedopt, inducedopt,roiopt)
% end
% return
%%
freq = 'evoked';
% save(save_name ,'Ybeta','Ytheta','Y','ltvall')
r = 0;
Yall = [];
YOut = [];
YCue = [];
YChoice = [];
ltv = [];
ltvOut = [];
ltvCue = [];
ltvChoice = [];
Ym = [];
YTheta = [];
YBeta = [];
YAlpha = [];
ntrials = [];
ntrialsOut = [];
ntrialsChoice = [];
ntrialsCue = [];
for sn = [1:9,11,12,14:16] %[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = cell(1);
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    
    for runs = 1:length(data_names)  
        data_name = data_names{runs};
        n = str2double(data_name(end-3));
        if isnan(n)
            n = 1;
        end
        
        load(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/',sub,'_',num2str(n)])
%         if isempty(ltv)
%             ltv = ltvall;
%         else
%             ltv(end+(1:size(ltvall,1)),:) = ltvall;
%         end
        % Evoked responses 200Hz sampling
%         Yall = cat(3,Yall,Y);
        % Induced responses 100Hz sampling
%         YTheta = cat(3,YTheta,Ytheta);
%         YBeta = cat(3,YBeta,Ybeta);
%         YAlpha = cat(3,YAlpha,Yalpha);
        % Check sign of evoked responses
%         Ym = cat(3,Ym,mean(Y,3));
%         ntrials = cat(1,ntrials,size(Y,3)); 

        if ~strcmp(freq, 'evoked')
            Ycue = eval(['Ycue_',freq]);
            Ychoice = eval(['Ychoice_',freq]);
            Yout = eval(['Yout_',freq]);
            ltvcue = eval(['ltvcue_',freq]);
            ltvchoice = eval(['ltvchoice_',freq]);
            ltvout = eval(['ltvout_',freq]);           
        end

        r = r+1;
        ltvout.recording = repmat(r,size(ltvout,1),1);
        ltvcue.recording =  repmat(r,size(ltvcue,1),1);
        ltvchoice.recording =  repmat(r,size(ltvchoice,1),1);
        
        YCue = cat(3,YCue,Ycue);
        YChoice = cat(3,YChoice,Ychoice);
        YOut = cat(3,YOut,Yout);
        % Check sign of evoked response
        Ym = cat(3,Ym, cat(2,mean(Ycue,3),mean(Yout,3),mean(Ychoice,3)));
%         Ym = cat(3,Ym, cat(2,mean(Yout,3),mean(Ychoice,3)));
        ntrialsChoice = cat(1,ntrialsChoice,size(Ychoice,3)); 
        ntrialsCue = cat(1,ntrialsCue,size(Ycue,3)); 
        ntrialsOut = cat(1,ntrialsOut,size(Yout,3)); 
        
        if isempty(ltvOut)
            ltvOut = ltvout;
            ltvCue = ltvcue;
            ltvChoice = ltvchoice;
        else
            ltvOut(end+(1:size(ltvout,1)),:) = ltvout;
            ltvCue(end+(1:size(ltvcue,1)),:) = ltvcue;
            ltvChoice(end+(1:size(ltvchoice,1)),:) = ltvchoice;
        end
        
        
    end
end

Ym0 = Ym;
% ltv = flipud(ltv);
clear ltvchoice ltvcue ltvout Y Yalpha Ybeta Ytheta Ycue Yout Ychoice
%% Sign flip for source localized evoked responses
% S = zeros(size(Yall));
if strcmp(freq,'evoked')
Sout = zeros(size(YOut));
Scue = zeros(size(YCue));
Schoice = zeros(size(YChoice));

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
    
%     time = linspace(-0.5,1,size(Yav,2));
%     figure; plot(time,mean(Yav,1))
  
%     zz = 0;
%     for z = 1:size(Yu,1)
%         S(ii,:,zz+(1:ntrials(z))) = s(z);
%         zz = zz + ntrials(z);
%     end
    
    zzo = 0;
    zzcu = 0;
    zzch = 0;
    
    for z = 1:size(Yu,1)
        Sout(ii,:,zzo+(1:ntrialsOut(z))) = s(z);
        zzo = zzo + ntrialsOut(z);
        Scue(ii,:,zzcu+(1:ntrialsCue(z))) = s(z);
        zzcu = zzcu + ntrialsCue(z);
        Schoice(ii,:,zzch+(1:ntrialsChoice(z))) = s(z);
        zzch = zzch + ntrialsChoice(z);
    end
    
    
end

% Yall = Yall.*S;


figure; set(gcf,'color','w')
subplot(231); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
subplot(232); imagesc(mean(YChoice,3)); caxis([-1 1]*4e-12)
subplot(233); imagesc(mean(YOut,3)); caxis([-1 1]*4e-12)

YCue = YCue.*Scue;
YChoice = YChoice.*Schoice;
YOut = YOut.*Sout;

subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
subplot(235); imagesc(mean(YChoice,3)); caxis([-1 1]*4e-12)
subplot(236); imagesc(mean(YOut,3)); caxis([-1 1]*4e-12)


% Flip ROIs to match each other
Yu = cat(2,mean(YCue,3),mean(YChoice,3),mean(YOut,3));
% Yu = cat(2,mean(YChoice(:,100:300,:),3),mean(YOut(:,80:200,:),3));

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
npoints = size(Scue,2);
Scue = repmat(s,[1,npoints,sum(ntrialsCue)]);
npoints = size(Schoice,2);
Schoice = repmat(s,[1,npoints,sum(ntrialsChoice)]);

YCue = YCue.*Scue;
YChoice = YChoice.*Schoice;
YOut = YOut.*Sout;

subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
subplot(235); imagesc(mean(YChoice,3));  caxis([-1 1]*4e-12)
subplot(236); imagesc(mean(YOut,3));  caxis([-1 1]*4e-12)
end


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

YCue = permute(YCue,[2,1,3]);
YCue = reshape(YCue,size(YCue,1)*size(YCue,2),size(YCue,3)); % nrois * npoints * ntrials

YChoice = permute(YChoice,[2,1,3]);
YChoice = reshape(YChoice,size(YChoice,1)*size(YChoice,2),size(YChoice,3)); % nrois * npoints * ntrials


%% Write data

dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_',freq,'_outcome.txt'],YOut);
dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_',freq,'_choice.txt'],YChoice);
dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_',freq,'_cue.txt'],YCue);

writetable(ltvOut,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_',freq,'_outcome.csv']);
writetable(ltvChoice,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_',freq,'_choice.csv']);
writetable(ltvCue,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_',freq,'_cue.csv']);
