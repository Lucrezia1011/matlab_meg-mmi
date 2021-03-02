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
Nlist([10,24]) = [];
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

% Define RPE groups at the group level 
R = cell(1,56);
Rlta = cell(1,56);
Nlistx = [1:9,11:56];  % only exclude subject 10
for sn = Nlistx
    sdan = num2str(meginfo.SDAN(sn));
    [R{sn},Rlta{sn}]= mmiBehavioralAll(sdan);
end
Rx = cell2mat(R');
Rltax = cell2mat(Rlta');
Nk = 3;
idx = kmeans([Rx,Rltax,abs(Rx)],Nk,'Replicates',30);

figure; hold all
for ii = 1:Nk
scatter(Rx(idx==ii),Rltax(idx==ii))
end
xlabel('RPE'); ylabel('RPE LTA')
grid on

s = zeros(1,Nk);
for ii = 1:Nk
    s(ii) = sum([mean(Rx(idx==ii)),mean(Rltax(idx==ii))]);
end
[~,idxs] = sort(s); % in increasing order
idx(idx==idxs(1)) = -1;
idx(idx==idxs(2)) = 0;
idx(idx==idxs(3)) = 1;

zz = 0;
for sn = Nlistx
    N = length(R{sn});
    R{sn} = idx(zz+(1:N));
    zz = zz+N;
end

filter_type = 'han'; % 'multi'; 'han'

for ii = 1:length(data_list)   
    sub = data_list{ii}(5:9);
    sn = (str2double(sub)==meginfo.SDAN);
%     mmiSensTFSprep(data_list{ii},filter_type,R{sn}) % run multitaper before
%     mmiSensTFSprep(data_list{ii},filter_type,[]) % get same number if wins and losses
end
return

%% Average over subjects
load /data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat

sn = 1;
sdan = num2str(meginfo.SDAN(sn));
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sdan,'/meg/'];
cd(data_path)
data_name = ['sub-',sdan,'_task-mmi3_run-1_meg.ds'];
processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/',data_name(1:end-3),'/'];

load(sprintf('%s/TFS_sens_filter-%s',processing_folder,filter_type));

[~,iia,iib] = intersect(channels,TFSwin.label);
% channels = TFSwin.label;
TFSwin_all = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time),length(Nlist) ]);
TFSlose_all =TFSwin_all;
TFSneut_all = TFSwin_all;

zz = 1;
base_sub = mean(TFScertain.powspctrm,3);
TFSwin_all(iia,:,:,zz) = TFSwin.powspctrm(iib,:,:) - base_sub(iib,:);
TFSlose_all(iia,:,:,zz) = TFSlose.powspctrm(iib,:,:) - base_sub(iib,:);
TFSneut_all(iia,:,:,zz) = TFSneut.powspctrm(iib,:,:) - base_sub(iib,:);

Nall = [TFSwin.ntrials TFSneut.ntrials TFSlose.ntrials];

for sn = Nlist(2:end) 
    
    sdan = num2str(meginfo.SDAN(sn));
    data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sdan,'/meg/'];
    cd(data_path)
    TFSwin_sub = zeros(size(TFSwin_all(:,:,:,1)));
    TFSlose_sub = zeros(size(TFSwin_sub));
    TFSneut_sub = zeros(size(TFSwin_sub));
    TFScertain_sub = zeros(size(TFSwin_sub));
    nWin = 0;
    nLose = 0;
    nCertain = 0;
    nNeut = 0;
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            
            processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/',data_name(1:end-3),'/'];           
            load(sprintf('%s/TFS_sens_filter-%s',processing_folder,filter_type));
                      
            if ~isempty(TFSwin)
                [~,iia,iib] = intersect(channels,TFSwin.label);
                TFSwin_sub(iia,:,:) = TFSwin_sub(iia,:,:) + TFSwin.powspctrm(iib,:,:)*TFSwin.ntrials;
                nWin = nWin + TFSwin.ntrials;
            end
            if ~isempty(TFSlose)
                [~,iia,iib] = intersect(channels,TFSlose.label);
                TFSlose_sub(iia,:,:) = TFSlose_sub(iia,:,:) + TFSlose.powspctrm(iib,:,:)*TFSlose.ntrials;
                nLose = nLose + TFSlose.ntrials;
            end
            if ~isempty(TFScertain)
                [~,iia,iib] = intersect(channels,TFScertain.label);
                TFScertain_sub(iia,:,:) = TFScertain_sub(iia,:,:) + TFScertain.powspctrm(iib,:,:)*TFScertain.ntrials;
                nCertain = nCertain + TFScertain.ntrials;
            end
            if ~isempty(TFSneut)
                [~,iia,iib] = intersect(channels,TFSneut.label);
                TFSneut_sub(iia,:,:) = TFSneut_sub(iia,:,:) + TFSneut.powspctrm(iib,:,:)*TFSneut.ntrials;
                nNeut = nNeut + TFSneut.ntrials;
            end
                        
        end
    end
    zz = zz +1;
    if nCertain ==0 || nWin ==0 || nLose ==0 || nNeut == 0
        warning('Sub %s: not enough trials!\n',sdan)
    else
        base_sub = mean(TFScertain_sub/nCertain,3);
        TFSwin_all(:,:,:,zz) = ((TFSwin_sub/nWin) - base_sub);
        TFSlose_all(:,:,:,zz) = ((TFSlose_sub/nLose) - base_sub);
        TFSneut_all(:,:,:,zz) = ((TFSneut_sub/nNeut) - base_sub);
    end
    Nall = Nall + [nWin nNeut nLose];
end

TFSwin.powspctrm = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSwin.label = channels;
TFSlose.powspctrm = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSlose.label = channels;
TFSneut.powspctrm =  zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSneut.label = channels;

% last subject has correct channels names
for ii = 1:length(channels)
    indC = abs(squeeze(TFSwin_all(ii,3,3,:)))>0; % ignore zeros and NaNs
    TFSwin.powspctrm(ii,:,:) = mean(TFSwin_all(ii,:,:,indC),4) ;
    
    indC = abs(squeeze(TFSlose_all(ii,3,3,:)))>0;
    TFSlose.powspctrm(ii,:,:) = mean(TFSlose_all(ii,:,:,indC),4) ;
    
    indC = abs(squeeze(TFSneut_all(ii,3,3,:)))>0;
    TFSneut.powspctrm(ii,:,:) = mean(TFSneut_all(ii,:,:,indC),4) ;
end

TFSwin.powspctrm = TFSwin.powspctrm*sqrt(Nall(1));
TFSneut.powspctrm = TFSneut.powspctrm*sqrt(Nall(2));
TFSlose.powspctrm = TFSlose.powspctrm*sqrt(Nall(3));

%% Average over trials
load /data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat

sn = 1;
sdan = num2str(meginfo.SDAN(sn));
data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sdan,'/meg/'];
cd(data_path)
data_name = ['sub-',sdan,'_task-mmi3_run-1_meg.ds'];
processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/',data_name(1:end-3),'/'];

load(sprintf('%s/TFS_sens_filter-%s',processing_folder,filter_type));

[~,iia,iib] = intersect(channels,TFSwin.label);
% channels = TFSwin.label;
TFSwin_all = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time),length(Nlist) ]);
TFSlose_all =TFSwin_all;
TFSneut_all = TFSwin_all;

zz = 1;
base_sub = mean(TFScertain.powspctrm,3);
base_sub = base_sub(iib,:);
base_sub = repmat(base_sub,[1,1,length(TFScertain.time)]);
TFSwin_all(iia,:,:,zz) = (TFSwin.powspctrm(iib,:,:) - base_sub)./ base_sub*(TFSwin.ntrials); % multiply the mean by sqrt(N) again to sum over all trials
TFSlose_all(iia,:,:,zz) = (TFSlose.powspctrm(iib,:,:)- base_sub)./ base_sub*(TFSlose.ntrials);
TFSneut_all(iia,:,:,zz) = (TFSneut.powspctrm(iib,:,:) - base_sub)./ base_sub*(TFSneut.ntrials);

Nall = [TFSwin.ntrials TFSneut.ntrials TFSlose.ntrials];

for sn = Nlist(2:end) 
    
    sdan = num2str(meginfo.SDAN(sn));
    data_path = ['/data/MBDU/MEG_MMI3/data/bids/sub-',sdan,'/meg/'];
    cd(data_path)
    TFSwin_sub = zeros(size(TFSwin_all(:,:,:,1)));
    TFSlose_sub = zeros(size(TFSwin_sub));
    TFSneut_sub = zeros(size(TFSwin_sub));
    TFScertain_sub = zeros(size(TFSwin_sub));
    nWin = 0;
    nLose = 0;
    nCertain = 0;
    nNeut = 0;
    
    for iiN = 1:3
        data_name = ['sub-',sdan,'_task-mmi3_run-',num2str(iiN),'_meg.ds'];
        if exist(data_name,'dir') && ~any(strcmp(data_name,data_exclude))
            
            processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sdan,'/',data_name(1:end-3),'/'];           
            load(sprintf('%s/TFS_sens_filter-%s',processing_folder,filter_type));
                      
            if ~isempty(TFSwin)
                [~,iia,iib] = intersect(channels,TFSwin.label);
                TFSwin_sub(iia,:,:) = TFSwin_sub(iia,:,:) + TFSwin.powspctrm(iib,:,:)*TFSwin.ntrials;
                nWin = nWin + TFSwin.ntrials;
            end
            if ~isempty(TFSlose)
                [~,iia,iib] = intersect(channels,TFSlose.label);
                TFSlose_sub(iia,:,:) = TFSlose_sub(iia,:,:) + TFSlose.powspctrm(iib,:,:)*TFSlose.ntrials;
                nLose = nLose + TFSlose.ntrials;
            end
            if ~isempty(TFScertain)
                [~,iia,iib] = intersect(channels,TFScertain.label);
                TFScertain_sub(iia,:,:) = TFScertain_sub(iia,:,:) + TFScertain.powspctrm(iib,:,:)*TFScertain.ntrials;
                nCertain = nCertain + TFScertain.ntrials;
            end
            if ~isempty(TFSneut)
                [~,iia,iib] = intersect(channels,TFSneut.label);
                TFSneut_sub(iia,:,:) = TFSneut_sub(iia,:,:) + TFSneut.powspctrm(iib,:,:)*TFSneut.ntrials;
                nNeut = nNeut + TFSneut.ntrials;
            end
                        
        end
    end
    zz = zz +1;
    base_sub = mean(TFScertain_sub/nCertain,3);
    base_sub = repmat(base_sub,[1,1,size(TFSwin_all,3)]);
    
    TFSwin_all(:,:,:,zz) = ((TFSwin_sub/nWin) - base_sub)./ base_sub*(nWin);
    TFSlose_all(:,:,:,zz) = ((TFSlose_sub/nLose) - base_sub)./ base_sub*(nLose);
    TFSneut_all(:,:,:,zz) = ((TFSneut_sub/nNeut) - base_sub)./ base_sub*(nNeut);
    
    Nall = Nall + [nWin nNeut nLose];
end

TFSwin.powspctrm = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSwin.label = channels;
TFSlose.powspctrm = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSlose.label = channels;
TFSneut.powspctrm = zeros([length(channels),length(TFSwin.freq),length(TFSwin.time)]);
TFSneut.label = channels;

% last subject has correct channels names
for ii = 1:length(channels)
    indC = abs(squeeze(TFSwin_all(ii,3,3,:)))>0; % ignore zeros and NaNs
    TFSwin.powspctrm(ii,:,:) = sum(TFSwin_all(ii,:,:,indC),4)/Nall(1) ;
    
    indC = abs(squeeze(TFSlose_all(ii,3,3,:)))>0;
    TFSlose.powspctrm(ii,:,:) = sum(TFSlose_all(ii,:,:,indC),4)/Nall(3);
    
    indC = abs(squeeze(TFSneut_all(ii,3,3,:)))>0;
    TFSneut.powspctrm(ii,:,:) = sum(TFSneut_all(ii,:,:,indC),4)/Nall(2) ;
end

TFSwin = rmfield(TFSwin,'ntrials');
TFSneut = rmfield(TFSneut,'ntrials');
TFSlose = rmfield(TFSlose,'ntrials');


%% Plot


cfg = [];
cfg.baseline     = 'no';%[TFScertain.time(1) TFScertain.time(end)];
cfg.baselinetype = 'absolute';
cfg.channel      = channels;
cfg.zlim         = [-1 1]*0.3;
cfg.xlim         = [-0.2 1];
cfg.showlabels   = 'yes';
cfg.layout       = 'CTF275_helmet.mat';
cfg.interactive  = 'no';

wind_size =  [491  66  1030   821];
figure(1); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSwin); 

figure(2); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSlose); 

figure(3); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSneut); 

TFSwin_neut = TFSwin;
TFSwin_neut.powspctrm = TFSwin.powspctrm - TFSneut.powspctrm;
figure(4); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSwin_neut); 

TFSlose_neut = TFSlose;
TFSlose_neut.powspctrm = TFSlose.powspctrm - TFSneut.powspctrm;
figure(5); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSlose_neut); 

TFSwin_lose = TFSlose;
TFSwin_lose.powspctrm = TFSwin.powspctrm - TFSlose.powspctrm;
figure(6); clf; set(gcf,'color','w','position', wind_size)
ft_multiplotTFR(cfg, TFSwin_lose); 

name_list = {'win';'lose';'neut';'win-neut';'lose-neut';'win-lose'};
for ii = 1:6
%     figure(ii);
%     saveas(gcf,sprintf('~/matlab/figures/TFS_trialav_%s.tif',name_list{ii}))
end
%%


cfg = [];
cfg.baseline     = 'no';%[TFScertain.time(1) TFScertain.time(end)];
cfg.baselinetype = 'absolute';
cfg.channel      = {'MRT21';'MRF35';'MRF25';'MRT11';'MRT32'};
cfg.zlim         = [-1 1]*0.3;
cfg.xlim         = [-0.2 1];
% cfg.showlabels   = 'yes';
% cfg.layout       = 'CTF275_helmet.mat';
cfg.interactive  = 'no';

figure; set(gcf,'color','w')
ft_singleplotTFR(cfg, TFSwin_neut); 
xlabel('Time (s)'); ylabel('Frequency (Hz)')