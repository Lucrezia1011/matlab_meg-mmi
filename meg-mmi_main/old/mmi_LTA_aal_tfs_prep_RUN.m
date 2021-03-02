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

% Rerun for hanning filter? 
% Check TFCE for positive and negative
roiopt = 'g'; filter_opt = 'han';
% for ii = 2:length(param_list) % redo 15
%     mmi_LTA_aal_tfs_prep(param_list{ii},roiopt,filter_opt)
% end
% return

% 
% base1 = mean(squeeze(mean(tfsout.powspctrm(:,ii,:,tfsout.time<0),1)),2);
% ii = 70;
% pcolor(tfsout.time,tfsout.freq,(squeeze(mean(tfsout.powspctrm(:,ii,:,:),1))-base1)./base1)
% shading flat
% colorbar
% colormap jet
% caxis([-.7 .7])

%%


YOut = [];
YCue = [];
YChoice = [];

ltvOut = [];
ltvCue = [];
ltvChoice = [];


r = 0;

for sn =[1:9,11,12,14:16] %[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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
        clear ltvchoice ltvcue ltvout Ycue Yout Ychoice 
        
        r = r+1;
        ltvout_tfs.recording = repmat(r,size(ltvout_tfs,1),1);
        ltvcue_tfs.recording =  repmat(r,size(ltvcue_tfs,1),1);
        ltvchoice_tfs.recording =  repmat(r,size(ltvchoice_tfs,1),1);
        
%         YCue = cat(1,YCue,Ycue.powspctrm);
%         YChoice = cat(1,YChoice,Ychoice.powspctrm);
%         YOut = cat(1,YOut,Yout.powspctrm);        
        
        base1 = mean(mean(tfscue.powspctrm(:,:,:,1:5),4),1);
        base2 = mean(mean(tfsout.powspctrm(:,:,:,1:5),4),1);
%         base1 = mean(mean(cat(1,Ycue.powspctrm,Ychoice.powspctrm,Yout.powspctrm),4),1);
%         base2 = base1; 
%         YCue = cat(1,YCue, (Ycue.powspctrm  - base1)./base1);
%         YChoice = cat(1,YChoice,(Ychoice.powspctrm - base1)./base1);
%         YOut = cat(1,YOut,(Yout.powspctrm - base2)./base2);        
        YCue = cat(1,YCue, tfscue.powspctrm );
        YChoice = cat(1,YChoice,tfschoice.powspctrm );
        YOut = cat(1,YOut,tfsout.powspctrm);   
        
        if isempty(ltvOut)
            ltvOut = ltvout_tfs;
            ltvCue = ltvcue_tfs;
            ltvChoice = ltvchoice_tfs;
        else
            ltvOut(end+(1:size(ltvout_tfs,1)),:) = ltvout_tfs;
            ltvCue(end+(1:size(ltvcue_tfs,1)),:) = ltvcue_tfs;
            ltvChoice(end+(1:size(ltvchoice_tfs,1)),:) = ltvchoice_tfs;
        end
        
        
    end
end

timeout = tfsout.time;
timecue = tfscue.time;
timechoice = tfschoice.time;
freqs = tfsout.freq;
clear ltvchoice_tfs ltvcue_tfs ltvout_tfs tfscue tfsout tfschoice

aal_labels = readcell('~/labels_AAL116_MNIv4.csv');

%%
% figure; set(gcf,'color','w'); 
clf
ii = 70; 
cax = 1e-23;
nneg = ltvCue.mood< (0.62 - 0.1);
npos = ltvCue.mood> (0.62 + 0.1);
% nneg = ltvCue.RPE == 0;
% npos = ltvCue.RPE > 4;
% nneg = ltvCue.RPE == 0;
% npos = ltvCue.RPE > 4;
set(gcf,'Name',aal_labels{ii})

subplot(321)
pcolor(timecue, freqs, squeeze(mean(YCue(npos,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Cue positive'))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])

subplot(322)
pcolor(timechoice, freqs, squeeze(mean(YChoice(npos,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Choice positive'))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])

subplot(323)
pcolor(timecue, freqs, squeeze(mean(YCue(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Cue negative '))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])

subplot(324)
pcolor(timechoice, freqs, squeeze(mean(YChoice(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Choice negative'))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])


subplot(325)
pcolor(timecue, freqs, squeeze(mean(YCue(npos,ii,:,:),1)) - squeeze(mean(YCue(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Cue in positive-negative'))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])

subplot(326)
pcolor(timechoice, freqs, squeeze(mean(YChoice(npos,ii,:,:),1)) -  squeeze(mean(YChoice(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Choice positive-negative'))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 40])
%%
% figure; set(gcf,'color','w'); 
clf
ii = 84;
cax  =1e-24;
% 
nneg = ltvOut.RPE<-2; % 116
npos = ltvOut.RPE>4 ; % 177
nneu = (abs(ltvOut.RPE)<1); % 219

% nneg = ltvOut.RPE_LTA<-4; % 159
% npos = ltvOut.RPE_LTA>7; % 153
% nneu = (abs(ltvOut.RPE_LTA)<1.8); % 152

subplot(321)
pcolor(timeout, freqs, squeeze(mean(YOut(npos,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Positive feedback in %s', aal_labels{ii}))
ylabel('frequency (Hz)'); ylim([1 50]);

subplot(323)
pcolor(timeout, freqs, squeeze(mean(YOut(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Negative feedback in %s', aal_labels{ii}))
ylabel('frequency (Hz)'); ylim([1 50]);

subplot(325)
pcolor(timeout, freqs, squeeze(mean(YOut(nneu,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Neutral feedback in %s', aal_labels{ii}))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 50]);
colormap jet


base1 = mean(squeeze(mean(YOut(nneu,ii,:,:),1)),2);
subplot(322)
pcolor(timeout, freqs, squeeze(mean(YOut(npos,ii,:,:),1)) - base1);
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Positive feedback in %s', aal_labels{ii}))
ylabel('frequency (Hz)'); ylim([1 50]);

subplot(324)
pcolor(timeout, freqs, squeeze(mean(YOut(nneg,ii,:,:),1))- base1);
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Negative feedback in %s', aal_labels{ii}))
ylabel('frequency (Hz)'); ylim([1 50]);

subplot(326)
pcolor(timeout, freqs, squeeze(mean(YOut(npos,ii,:,:),1))-squeeze(mean(YOut(nneg,ii,:,:),1)));
shading interp; colorbar; caxis([-1 1]*cax)
title(sprintf('Pos-Neg feedback in %s', aal_labels{ii}))
xlabel('time(s)'); ylabel('frequency (Hz)'); ylim([1 50]);
colormap jet

% 
%%
% Time x frequency X ROI x trial
YOut = permute(YOut,[4,3,2,1]);
YOut = reshape(YOut,size(YOut,1)*size(YOut,2)*size(YOut,3),size(YOut,4));

YCue = permute(YCue,[4,3,2,1]);
YCue = reshape(YCue,size(YCue,1)*size(YCue,2)*size(YCue,3),size(YCue,4)); 

YChoice = permute(YChoice,[4,3,2,1]);
YChoice = reshape(YChoice,size(YChoice,1)*size(YChoice,2)*size(YChoice,3),size(YChoice,4)); 

% get rid of some zeros
YOut = YOut*1e24;
YCue = YCue*1e24;
YChoice = YChoice*1e24;
%% Write data

dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_tfs_outcome.txt',YOut);
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_tfs_choice.txt',YChoice);
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/meg_trials_tfs_cue.txt',YCue);

writetable(ltvOut,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_tfs_outcome.csv']);
writetable(ltvChoice,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_tfs_choice.csv']);
writetable(ltvCue,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/latent_vars_tfs_cue.csv']);
