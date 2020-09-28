clear all
% close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')

%% Extract ECG component from all MMI data sets: compare with mood timecourse


subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

param_list = [];

zz= 0;
for sn = [1:9,11:12,14:16] % all subjects with continuos recordings and latent variables
    
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = cell(1);
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ...
                ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') && ...
                exist([data_path,name_list(ii).name,'/beamforming/mcg.mat/'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    if jj>0
        for runs = 1:length(data_names)
            zz = zz +1;
            param_list{zz} = data_names{runs};
        end
    end
end

%% Missing 
% close all; clc;
% for zz = 1:length(param_list)  % 10(bad quality), 16,18 missing
%     data_name = param_list{zz};
%     sub = data_name(1:5);
%     data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
%     cd(data_path)
%     processing_folder = [data_path,data_name,'/beamforming'];
%     
%     if ~exist([processing_folder,'/mcg.mat'],'file')
%     
%     load([processing_folder,'/ICA_artifacts.mat']);
%     f = 1200;
%     ecgm = ft_preproc_bandpassfilter(comps.trial{1}, comps.fsample, [1 35],[],'but');
%     l = comps.sampleinfo(2)/f;
%     figure; plot(ecgm'); xlim([0 f*10])
%     
%     ic = [];
%     for ii = 1:size(ecgm,1)
%         [pksm,~] = findpeaks(ecgm(ii,:),f,'MinPeakDistance',0.4,'MinPeakWidth',0.01,'MinPeakProminence',5e-12);
%         if size(pksm,2) > l/2
%             ic = cat(1,ic,ii);
%         end
%     end
%     
%     if ~isempty(ic)
%     [m,ici] = max(std(ecgm(ic,:)'));
%  
%     ic = ic(ici);
%     ecg = ecgm(ic,:);
%     if nnz(ecg > (mean(ecg)+2*m))  <  nnz(ecg < (mean(ecg)-2*m)) 
%         ecg = -ecg;
%     end
%     
%     [pks,locs] = findpeaks(ecg,f,'MinPeakDistance',0.4,'MinPeakWidth',0.01,'MinPeakProminence',4*m);%,'MaxPeakWidth',0.3);
% 
%     
%     figure;  set(gcf,'color','w')
%     subplot(211);
%     plot(linspace(1/f,l,l*f),ecg,'k')
%     hold on
%     plot(locs,pks,'*r')
%     ylabel('Magnetic field (T)')
%     title('MEG IC peaks')
%     xlabel('time(s)')
%     subplot(212);
%     plot((locs(1:end-1)+locs(2:end))/2, 60./diff(locs),'k')
%     xlabel('time(s)');
%     ylabel('Beats per minute (BPM)')
%     title('Heart rate (RR interval)')
%     
%     save([processing_folder,'/mcg.mat'],'locs')
%     end
%     end
% end  

    %% Read events
close all

c_mRR = zeros(1,14);
c_SDRR = zeros(1,14);
c_RMSSD = zeros(1,14);
c_PVLF = zeros(1,14);
c_PLF = zeros(1,14);
c_PHF = zeros(1,14);
c_TP = zeros(1,14);
c_PLHF = zeros(1,14);
c_SD1 = zeros(1,14);
c_SD2 = zeros(1,14);
c_SD12 = zeros(1,14);
s_mRR = zeros(1,14);
s_SDRR = zeros(1,14);
s_RMSSD = zeros(1,14);
s_PVLF = zeros(1,14);
s_PLF = zeros(1,14);
s_PHF = zeros(1,14);
s_TP = zeros(1,14);
s_moodvar = zeros(1,14);
s_moodbase = zeros(1,14);
SD1 = zeros(1,14);
SD2 = zeros(1,14);

zz= 0;
for sn = [1:9,11:12,14:16] % all subjects with continuos recordings and latent variables
    
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    name_list = dir;
    data_names = [];
    jj = 0;
    for ii = 1:length(name_list)
        if strncmp(name_list(ii).name, data_name, 18) && ...
                ~strcmp(name_list(ii).name, '24201MMI_mmi3_proc1.ds') && ...
                exist([data_path,name_list(ii).name,'/beamforming/mcg.mat'],'file')
            jj = jj+1;
            data_names{jj} = name_list(ii).name;
        end
    end
    zz = zz+1;
    
    Mood_1 = [];
    mRR_1 = [];
    SDRR_1 = [];
    RMSSD_1 = [];
    PVLF_1 = [];
    PLF_1 = [];
    PHF_1 = [];
    TP_1 = [];
    SD_1 = [];
    SD_2 = [];
    RR_1 = [];
    Mooda_1 = [];
    pocar = [];
    for runs = 1:length(data_names)
        
        data_name = data_names{runs};

        processing_folder = [data_path,data_name,'/beamforming'];
        load([processing_folder,'/mcg.mat']);
        f = 1200;
        bv_match = match_triggers_fc(data_name);
        
        mood_match = bv_match.ratemood;
        blockmood_match = bv_match.blockmood;
        tasktime = bv_match.time;
        
        %% Mood
        
        tasktime = tasktime - tasktime(1) + 1/f;
        
        indbm = blockmood_match.sample~=0;
        indm = mood_match.sample~=0;
        
        [x,ind] = sort([blockmood_match.sample(indbm), mood_match.sample(indm)]);
        v = [blockmood_match.mood(indbm), mood_match.mood(indm)];
        Fmood = griddedInterpolant(tasktime(x),v(ind),'pchip');
        
        % 30s window around mood report?
        
        %%
        
        time_mcg = (locs(1:end-1)+locs(2:end))/2;
        RR = diff(locs)*1e3;
        
        Fsamp = 2; % resample with a uniform 2Hz frequency
        
        F = griddedInterpolant(time_mcg,RR,'pchip');
        t1 = ceil(time_mcg(1)*Fsamp)/Fsamp;
        t2 = floor(time_mcg(end)*Fsamp)/Fsamp;
        t = t1:1/Fsamp:t2;
        RRi = F(t);
          
        
        mRR = zeros(1,length(x));
        SDRR = zeros(1,length(x));
        RMSSD = zeros(1,length(x));
        PVLF = zeros(1,length(x));
        PLF = zeros(1,length(x));
        PHF = zeros(1,length(x));
        TP = zeros(1,length(x));
        mood = v(ind);
        k_ind = [];
        
        
        % 45degrees rotation matrix
        R = [cos(pi/4) -sin(pi/4)
            sin(pi/4)  cos(pi/4)];
        sd1 = zeros(1,length(x));
        sd2 = zeros(1,length(x));
        
        for k = 1:length(x)
            
            w = 30/2;
            xx = x(k);
            timeind = round(tasktime(xx) - w):round(tasktime(xx) + w);
            [~,tt1] = min( abs( t - (tasktime(xx) - w)));
            [~,tt2] = min( abs( t - (tasktime(xx) + w)));
            

            RRk = RRi(tt1:tt2);
            if length(RRk) == w*Fsamp*2+1
                mRR(k) = mean(RRk);
                nfft = 256;
%                 Pk = abs(fft(RRk-mean(RRk),nfft));
                Pk = abs(fft(RRk-mean(RRk),nfft));
                Pk = Pk(1:nfft/2+1);
                freqsk = linspace(0,Fsamp/2,nfft/2+1);
%                 plot(freqsk,Pk)
%                 [Pk,freqsk] = pwelch(RRk-mean(RRk),10,[],256,Fsamp);
%                 [Pk,freqsk] = pwelch(RRk,2*Fsamp,[],256,Fsamp);
                % RMSDD root mean square successive differences
                RMSSD(k) = sqrt(mean(diff(RRk).^2));
                SDRR(k) = std(RRk);

                [~, ind1] = min(abs(freqsk - 0.04 ));
                [~, ind2] = min(abs(freqsk - 0.15 ));

                PLF(k) = trapz(freqsk(ind1:ind2),Pk(ind1:ind2) );

                [~, ind1] = min(abs(freqsk - 0.15 ));
                [~, ind2] = min(abs(freqsk - 0.40 ));

                PHF(k) = trapz(freqsk(ind1:ind2),Pk(ind1:ind2) );            
                TP(k) = trapz(freqsk,Pk);
                
                % Should use the non interpolated R-R?
                pocart = R * [RRk(1:end-1);RRk(2:end)];
        %         figure; scatter(pocart(1,:),pocart(2,:))
                
                sd1(k) = std(pocart(1,:));
                sd2(k) = std(pocart(2,:));

                k_ind = [k_ind,k];
            end
        end
          
        Mood_1 = cat(2,Mood_1,mood(k_ind));
        mRR_1 = cat(2,mRR_1,mRR(k_ind));
        SDRR_1 = cat(2,SDRR_1,SDRR(k_ind));
        RMSSD_1 = cat(2,RMSSD_1,RMSSD(k_ind));
        
        PLF_1 = cat(2,PLF_1,PLF(k_ind));
        PHF_1 = cat(2,PHF_1,PHF(k_ind));
        TP_1 = cat(2,TP_1,TP(k_ind));
        SD_1 = cat(2,SD_1,sd1(k_ind));
        SD_2 = cat(2,SD_2,sd2(k_ind));

        Fsamp = 1; % resample with a uniform 1Hz frequency
        
        F = griddedInterpolant(time_mcg,RR,'pchip');
        t1 = ceil(time_mcg(1)*Fsamp)/Fsamp;
        t2 = floor(time_mcg(end)*Fsamp)/Fsamp;
        t = t1:1/Fsamp:t2;
        RRi = F(t);
        RR_1 = cat(2,RR_1,RRi);
        pocar = cat(2,[RRi(1:end-1);RRi(2:end)]);

        % Full mood timecourse 
        Mooda_1 = cat(2,Mooda_1,mood);
    end
    
    % figure; set(gcf,'color','w')
    % yyaxis left; plot(tasktime(x),[mRR;SDRR;RMSSD]); ylabel('ms')
    % yyaxis right; plot(tasktime(x),[PVLF;PLF;PHF;TP]); ylabel('ms^2')
    % xlabel('time (s)')
    % legend('Mean RR (ms)','STDRR (ms)','RMSSD (ms)','Power VLF ms^2',...
    %     'Power LF ms^2','Power HF ms^2','Total power ms^2')
    corrtype = 'pearson';
    c_mRR(zz) = corr(Mood_1',mRR_1','type',corrtype);
    c_SDRR(zz) = corr(Mood_1',SDRR_1','type',corrtype);
    c_RMSSD(zz) = corr(Mood_1',RMSSD_1','type',corrtype);
    
    c_PLF(zz) = corr(Mood_1',PLF_1','type',corrtype);
    c_PHF(zz) = corr(Mood_1',PHF_1','type',corrtype);
    c_TP(zz) = corr(Mood_1',TP_1','type',corrtype);
    c_PLHF(zz) = corr(Mood_1',(PLF_1./PHF_1)','type',corrtype);
    
    c_SD1(zz) = corr(Mood_1',SD_1','type',corrtype);
    c_SD2(zz) = corr(Mood_1',SD_2','type',corrtype);
    c_SD12(zz) = corr(Mood_1',(SD_1./SD_2)','type',corrtype);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
               
        nfft = 1024;
        % Low frequency power% correlates with mood only when baseline
        % correcting R-R before frequency analysis
        [P,freqs] = pwelch(RR_1-mean(RR_1),256*Fsamp,[],nfft,Fsamp);
%         [P,freqs] = pwelch(RR_1,256*Fsamp,[],nfft,Fsamp);

        % figure; set(gcf,'position', [978     9   559   906])
        %
        % subplot(411)
        % plot(t, RRi)
        % hold on
        % xlabel('time(s)'); ylabel('R-R interval (ms)')
        % title('Heart rate (R-R interval)')
        %
        %
        % subplot(412)
        % plot(freqs, P)
        %
        % xlabel('Frequency (Hz)')
        % ylabel('Power (ms^2/Hz)')
        % set(gcf,'color','w')
        % xlim([0 0.5])
        % title('Heart rate PSD')
     
        % RMSDD root mean square successive differences
        s_RMSSD(zz) = sqrt(mean(diff(RR_1).^2));
        s_SDRR(zz) = std(RR_1);
       
        [~, ind1] = min(abs(freqs - 0.0033 ));
        [~, ind2] = min(abs(freqs - 0.04 ));
        
        s_PVLF(zz) = trapz(freqs(ind1:ind2),P(ind1:ind2) );
   
        [~, ind1] = min(abs(freqs - 0.04 ));
        [~, ind2] = min(abs(freqs - 0.15 ));
        
        s_PLF(zz) = trapz(freqs(ind1:ind2),P(ind1:ind2) );
        
        [~, ind1] = min(abs(freqs - 0.15 ));
        [~, ind2] = min(abs(freqs - 0.40 ));
        
        s_PHF(zz) = trapz(freqs(ind1:ind2),P(ind1:ind2) );
        
        s_TP(zz) = trapz(freqs,P);
        s_mRR(zz) = mean(RR_1);        

        s_moodvar(zz) = var(Mooda_1);
        s_moodbase(zz) = Mooda_1(1);
%         clc
%         fprintf(['Mean RR %.1fms, STDRR %.1fms and RMSSD %.1fms\n',...
%             'Power VLF %.2f ms^2, LF %.2f ms^2, HF %.2f ms^2\nTotal power %.2f ms^2\n\n'],...
%             mean(RR_1),SDRRs,RMSSDs,PVLFs,PLFs,PHFs,TPs)
        
        %
       
%       Poincare measures
    
%         figure
%         scatter(pocar(1,:),pocar(2,:))
%         hold on; plot(650:1000,650:1000)
%         xlabel('R-R_t (ms)'); ylabel('R-R_{t+1} (ms)')
%         title('Poincare plot of RR intervals')
%         axis equal; %([750 1500 750 1500])

        % 45degrees rotation matrix
        R = [cos(pi/4) -sin(pi/4)
          sin(pi/4)  cos(pi/4)];

        pocart = R * pocar;
%         figure; scatter(pocart(1,:),pocart(2,:))

        SD1(zz) = std(pocart(1,:));
        SD2(zz) = std(pocart(2,:));

end

%%
c_all = [c_mRR;c_SDRR;c_RMSSD;c_PLF;c_PHF;c_TP;c_PLHF;c_SD1;c_SD2;c_SD12];
figure; set(gcf,'color','w')
errorbar(1:size(c_all,1),mean(c_all,2),std(c_all,0,2),'o')
xlim([0 size(c_all,1)+1])
set(gca,'XTick',1:size(c_all,1),'XTickLabels',...
    {'mRR';'SDRR';'RMSSD';'PLF';'PHF';'TP';'PLF/PHF';'SD1';'SD2';'SD1/SD2'})
ylabel('Pearson correlation (r)')
title('Correlation of Heart rate variability parameters with mood')

mdds = [1,4,5,13,14];
hvs = [2,3,6:12];
sall = [s_mRR; s_SDRR; s_RMSSD; s_PVLF; s_PVLF./s_TP; s_PLF; (s_PLF./s_TP); ...
    s_PHF; (s_PHF./s_TP); (s_PLF./s_PHF); SD1; SD2; (SD1./SD2)];
T = table;
T.HV_mean = mean(sall(:,hvs),2);
T.HV_stdev = std(sall(:,hvs),0,2);
T.MDD_mean = mean(sall(:,mdds),2);
T.MDD_stdev = std(sall(:,mdds),0,2);

T.Properties.RowNames = {'R-R (ms)'; 'SDRR (ms)'; 'RMSSD (ms)'; 'PVLF (ms^2)'; ...
    'PVLF (%)'; 'PLF (ms^2)'; 'PLF (%)'; ...
    'PHF (ms^2)'; 'PHF (%)'; 'PLF/PHF'; 'SD1 (ms)'; 'SD2 (ms)'; 'SD1/SD2'};



corrtype = 'spearman';
c = zeros(2,10);
c(2,1) = corr(s_moodvar',s_SDRR','type',corrtype);
c(2,2) = corr(s_moodvar',s_RMSSD','type',corrtype);
c(2,3) = corr(s_moodvar',s_PLF','type',corrtype);
c(2,4) = corr(s_moodvar',(s_PLF./s_TP)','type',corrtype);
c(2,5) = corr(s_moodvar',s_PHF','type',corrtype);
c(2,6) = corr(s_moodvar',(s_PHF./s_TP)','type',corrtype);
c(2,7) = corr(s_moodvar',(s_PLF./s_PHF)','type',corrtype);
c(2,8) = corr(s_moodvar',SD1','type',corrtype);
c(2,9) = corr(s_moodvar',SD2','type',corrtype);
c(2,10) = corr(s_moodvar',(SD1./SD2)','type',corrtype);

c(1,1) = corr(s_moodbase',s_SDRR','type',corrtype);
c(1,2) = corr(s_moodbase',s_RMSSD','type',corrtype);
c(1,3) = corr(s_moodbase',s_PLF','type',corrtype);
c(1,4) = corr(s_moodbase',(s_PLF./s_TP)','type',corrtype);
c(1,5) = corr(s_moodbase',s_PHF','type',corrtype);
c(1,6) = corr(s_moodbase',(s_PHF./s_TP)','type',corrtype);
c(1,7) = corr(s_moodbase',(s_PLF./s_PHF)','type',corrtype);
c(1,8) = corr(s_moodbase',SD1','type',corrtype);
c(1,9) = corr(s_moodbase',SD2','type',corrtype);
c(1,10) = corr(s_moodbase',(SD1./SD2)','type',corrtype);

figure;  set(gcf,'color','w')
bar(c); grid on; ylim([-0.8 0.8])
legend('SDRR','RMSSD','PLF','PLF%','PHF','PHF%','PLF/PHF','SD1','SD2','SD1/SD2',...
    'Location','bestoutside')
ylabel('Spearman correlation (r)')
title('Correlation of Heart rate variability parameters with mood')
set(gca,'XTick',1:2,'XTicklabels',{'Baseline mood';'Mood variance'})


return
%%

cfg.dataset     = 'MEG_ecg2.ds';
cfg.bpfilter = 'yes';
cfg.bpfreq  =[1 35];
ecg = ft_preprocessing(cfg);

f = ecg.fsample;
l = ecg.sampleinfo(2)/f; % recording length (s)

figure; set(gcf,'color','w')
subplot(211)
plot(linspace(1/f,l,l*f),(ecg_raw.trial{1}),'k')
ylabel('Electic potential (V)')
title('ECG signal 24/1 highpass filter 0.5Hz')

subplot(212)
plot(linspace(1/f,l,l*f),(ecg.trial{1}),'k')
xlabel('time(s)')
ylabel('Electic potential (V)')
title('ECG signal 31/1 bandpass filter 1-35Hz')

dlmwrite('/data/liuzzil2/ecg_day1.txt',[0:(1/ecg.fsample):300-(1/ecg.fsample);ecg_raw.trial{1}*1e3]','delimiter','\t','precision','%3.6f');
dlmwrite('/data/liuzzil2/ecg_day2.txt',[0:(1/ecg.fsample):300-(1/ecg.fsample);ecg.trial{1}*1e3]','delimiter','\t','precision','%3.6f');

%%

Fsamp = 1; % resample with a uniform 1Hz frequency

% 10.5-13 Hz high alpha,  13-20Hz low beta
data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [4 8],[],'but');
ht = hilbert(data_filt');
env = abs(ht);

envd = resample(env,Fsamp,data_clean.fsample);

F = griddedInterpolant(time_mcg,RR,'pchip');

t1 = ceil(time_mcg(1)*Fsamp)/Fsamp;
t2 = floor(time_mcg(end)*Fsamp)/Fsamp;
t = t1:1/Fsamp:t2;
RRi = F(t);

% ULF 0.0033-0.04,  LF 0.04-0.15, HF 0.15-0.4
RRULF = ft_preproc_bandpassfilter(RRi, Fsamp, [0.0033 0.04],[],'but');
RRLF = ft_preproc_bandpassfilter(RRi, Fsamp, [0.04 0.15],[],'but');
RRHF = ft_preproc_bandpassfilter(RRi, Fsamp, [0.15 0.4],[],'but');

% RRULF = abs(hilbert(RRULF));
RRLF = abs(hilbert(RRLF));
RRHF = abs(hilbert(RRHF));
RRULF = abs(hilbert(RRULF));


figure; plot(t,RRHF)
hold on
plot(t,RRLF)
plot(t,RRULF)


envd = envd(Fsamp*t1:(Fsamp*t2),:);

C = corr(zscore(RRi)',zscore(envd));
Cu = corr(zscore(RRULF)',zscore(envd));
Cl = corr(zscore(RRLF)',zscore(envd));
Ch = corr(zscore(RRHF)',zscore(envd));


figure;
subplot(131)
cz = data;
cz.time{1} = 1;
cz.trial{1} = Cu';
cz.avg = Cu';
cz.sampleinfo = [1 1];
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.zlim = [-.3 .3];
ft_topoplotER(cfg,cz)
title('Corr of ULF of RR and \alpha_2')

subplot(132)
cz.trial{1} = Cl';
cz.avg = Cl';
ft_topoplotER(cfg,cz)
title('Corr of LF of RR and \alpha_2')

subplot(133)
cz.trial{1} = Ch';
cz.avg = Ch';
ft_topoplotER(cfg,cz)
title('Corr of HF of RR and \alpha_2')



%% Try on MMI data

sub = '24213';
cd(['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub])
load([sub,'MMI_mmi3_proc.ds/beamforming/ICA_artifacts.mat'])

figure; plot(comps.trial{1}')

ecgm = ft_preproc_bandpassfilter(-comps.trial{1}(1,:), comps.fsample, [1 35],[],'but');
[pksm,locsm] = findpeaks(ecgm,f,'MinPeakDistance',0.5,'MinPeakWidth',0.01,'MinPeakProminence',2e-12);
l = comps.sampleinfo(2)/f;
figure;  set(gcf,'color','w')
subplot(311)
plot(linspace(1/f,l,l*f),ecgm,'k')
hold on
plot(locsm,pksm,'*r')
ylabel('Magnetic field (T)')
title([sub,' MEG IC peaks'])
xlabel('time(s)')
xlim([0 100])


Fsamp = 1; % resample with a uniform 10Hz frequency

RR = diff(locsm)*1e3;

locsdm = (locsm(1:end-1)+locsm(2:end))/2;

F = griddedInterpolant(locsdm,RR,'pchip');
t1 = ceil(locsdm(1)*Fsamp)/Fsamp;
t2 = floor(locsdm(end)*Fsamp)/Fsamp;
t = t1:1/Fsamp:t2;
RRi = F(t);

% Fbps = fft(bps - mean(bps),2048*2); % mean correct to eliminate 
% P = abs(Fbps(1:ceil(length(Fbps)/2)));
% freqs= linspace(0,Fsamp/2,length(P));
nfft = 1024;
[P,freqs] = pwelch(RRi-mean(RRi),256*Fsamp,[],nfft,Fsamp);

subplot(312)
plot(t, RRi, 'k')
xlabel('time(s)');
ylabel('R-R interval (ms)')
title('Heart rate (RR interval)')
axis tight
% ylim([0.6 1])
subplot(313)
plot(freqs, P ,'k')

xlabel('Frequency (Hz)')
ylabel('Power (ms^s/Hz)')
xlim([0 0.5]); ylim([0 1]*5e4)
title('Heart rate PSD')


% RMSDD root mean square successive differences
RMSSD = sqrt(mean(diff(RR).^2));
SDRR = std(RR);

[~, ind1] = min(abs(freqs - 0.0033 ));
[~, ind2] = min(abs(freqs - 0.04 ));

PVLF = trapz(freqs(ind1:ind2),P(ind1:ind2) );

[~, ind1] = min(abs(freqs - 0.04 ));
[~, ind2] = min(abs(freqs - 0.15 ));

PLF = trapz(freqs(ind1:ind2),P(ind1:ind2) );

[~, ind1] = min(abs(freqs - 0.15 ));
[~, ind2] = min(abs(freqs - 0.40 ));

PHF = trapz(freqs(ind1:ind2),P(ind1:ind2) );

TP = trapz(freqs,P );

clc
fprintf([sub,':\nMean RR %.1fms, STDRR %.1fms and RMSSD %.1fms\n',...
    'Power VLF %.2f ms^2, LF %.2f ms^2, HF %.2f ms^2\nTotal power %.2f ms^2\n\n'],...
    mean(RR),SDRR,RMSSD,PVLF,PLF,PHF,TP)


figure; set(gcf,'color','w')
scatter(RR(1:end-1),RR(2:end))
xlabel('R-R_t (ms)'); ylabel('R-R_{t+1} (ms)')
title([sub,': Poincare plot of RR intervals'])
axis equal
