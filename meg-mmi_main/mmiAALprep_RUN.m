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


analy_case = 'confirm';
roiopt = 'sens'; % running for grid
switch roiopt
    case 'grid'
        subexclude = [subexclude,26,49,53];
end

switch analy_case
    case 'confirm'
        subexclude = [subexclude,1:12,14:16];
    case 'explore'
        subexclude = [subexclude,13,17:56];
end

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

% for ii = 1:length(data_list) 
%     data_name = data_list{ii};
%     sub = data_name(5:9);
%     processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
%     mmiAALprep(data_list{ii},twind,roiopt)
% 
% end
% return
%%
freq = 'outcome';% 'outcome', 'cue' , 'choice'
% save(save_name ,'Ybeta','Ytheta','Y','ltvall')
r = 0;
Yall = [];
ltv = [];

Ym = [];

ntrials = [];
ntrialsOut = [];

if strcmp(roiopt,'sens')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    cd('/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/')
    hdr = ft_read_header(data_list{1});
    channelsall = hdr.label(strcmp(hdr.chantype,'meggrad'));
end
time = linspace(twind(1),twind(2),360);

for sn = 1:length(data_list) % all subjects with continuos recordings and latent variables
    
    data_name = data_list{sn};   
    sub = data_name(5:9);
    processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save_name = [processing_folder,'evoked_',freq,'_',roiopt,'.mat'];
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
        
        
        T = struct;
        T.label = channelSub;
        T.sampleinfo = [1 size(Yout,2)];
        T.chantype(1:length(T.label),1) = {'meggrad'    };
        T.grad = hdr.grad;

        T.trial = cell(1,ntrials);
        for n = 1:ntrials
            T.trial{n} = Yout(:,:,n);
        end
        T.time(1:ntrials) = {time};

        cfg =[];
        cfg.neighbours = neighbours;
        cfg.channel = T.label;
        T1 = ft_megplanar(cfg,T);

        cfg =[];
        cfg.demean = 'yes';
        T2 = ft_combineplanar([],T1);    
        
        % plot to check
        [timelock] = ft_timelockanalysis([], T);
        [timelockp] = ft_timelockanalysis([], T2);
        
        subplot(1,2,1)
        cfg = [];
        cfg.layout = 'CTF275_helmet.mat';
        cfg.xlim = [.1 .2];
        cfg.zlim = [-1 1]*5e-14;
        ft_topoplotER(cfg, timelock);
        
        subplot(1,2,2)
        cfg.zlim = [0 2]*1e-13;
        ft_topoplotER(cfg, timelockp)
        
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
    
    
    
end

Ym0 = Ym;
% ltv = flipud(ltv);
clear ltvchoice ltvcue ltvout Y Yalpha Ybeta Ytheta Ycue Yout Ychoice
%% Sign flip for source localized evoked responses
if strcmp(roiopt,'AAL')
% S = zeros(size(Yall));

Sout = zeros(size(Yall));
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
subplot(131); imagesc(mean(Yall,3)); caxis([-1 1]*.3)
title(sprintf('Grandaverage over all subjects\n of beamformed evoked responses'))
xlabel('timepoint'); ylabel('ROI')

% YCue = YCue.*Scue;
% YChoice = YChoice.*Schoice;
Yall = Yall.*Sout;

% subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
% subplot(235); imagesc(mean(YChoice,3)); caxis([-1 1]*4e-12)
subplot(132); imagesc(mean(Yall,3)); caxis([-1 1]*.3)
title(sprintf('Maximised ROI temporal\n correlation over subject'))
xlabel('timepoint'); ylabel('ROI')
% Flip ROIs to match each other
% Yu = cat(2,mean(YCue,3),mean(YChoice,3),mean(YOut,3));
Yu = mean(Yall,3);

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
Yall = Yall.*Sout;

% subplot(234); imagesc(mean(YCue,3)); caxis([-1 1]*4e-12)
% subplot(235); imagesc(mean(YChoice,3));  caxis([-1 1]*4e-12)
subplot(133); imagesc(mean(Yall,3));  caxis([-1 1]*.3)
title(sprintf('Maximised temporal correlation\n over all ROIs'))
xlabel('timepoint'); ylabel('ROI'); 
c = colorbar;
c.Position = [0.93 c.Position(2:4)];
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

%% Write data

switch roiopt
    case 'AAL'
        if strcmp(freq,'outcome')
            
            Yall = permute(Yall,[2,1,3]);
            Yall = reshape(Yall,size(Yall,1)*size(Yall,2),size(Yall,3)); % nrois * npoints * ntrials

            dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/meg_trials_evoked_',freq,'.txt'],Yall);
            writetable(ltv,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/latent_vars_evoked_',freq,'.csv']);
        else
            subs = unique(ltv.subject);
            Ym = zeros(size(Yall,1),size(Yall,2),length(subs));
            for s = 1:length(subs)
                Ym(:,:,s) = mean(Yall(:,:,ltv.subject == subs(s)),3);
            end
            Ym = permute(Ym,[2,1,3]);
            Ym = reshape(Ym,size(Ym,1)*size(Ym,2),size(Ym,3)); % nrois * npoints * ntrials

            dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/meg_trials_evoked_',freq,'.txt'],Ym);
        end
    case 'sens'
        writetable(ltv,['/data/MBDU/MEG_MMI3/results/mmiTrial_sens/evoked_outcome/latent_vars_evoked_',freq,'.csv']);
        
        Yall = permute(Yall,[2,3,1]);
        nrois =  size(channels,1);  
        
        for nn = 1:nrois
            n = num2str(nn);
            if size(n,2) == 1
                n = ['00',n];
            elseif size(n,2) == 2
                n = ['0',n];
            end
            dlmwrite(sprintf(['/data/MBDU/MEG_MMI3/results/mmiTrial_sens/',...
                'evoked_outcome/meg_trials/sens_%s.txt'],n),Yall(:,:,nn));
            
        end
end

%% Save averages
dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_sens/evoked_time/meg_trials_evoked_',freq,'.txt'],mean(Yall,3));

time = linspace(-.5,.7,360); % linspace(-.2,1,360)
M=mean(Yall(:,:,ltv.E>median(ltv.E)),3) - mean(Yall(:,:,ltv.E<median(ltv.E)),3);
v = squeeze(var(M,0,2));
v = max(abs(M),[],2);
[~,ind]=sort(v,'descend');
figure;clf; set(gcf,'color','w')
for c = 1:2
chans = ind(c);

subplot(2,1,c)
plot(time,mean(Yall(chans,:,ltv.E<median(ltv.E)),3)')
hold on
plot(time,mean(Yall(chans,:,ltv.E>median(ltv.E)),3)')
hold on
plot(time,mean(Yall(chans,:,ltv.E>median(ltv.E)),3)'-...
    mean(Yall(chans,:,ltv.E<median(ltv.E)),3)','k')
grid on; ylim([-1 1]*.5e-13)
title(['sensor ',channels{chans}])
end
legend('Low expectation','high expectation','diff')

v = squeeze(var(Ym,0,2));
[~,ind]=sort(v,'descend');
u=unique(ind(1:10,:));
N = zeros(1,length(u));
for iiU = 1:length(u)
   N(iiU) = nnz(ind(1:10,:)==u(iiU) );
end
[~,indN] = sort(N,'descend');
chans = sort(u(indN(1:10)));

figure
subplot(211)
plot(linspace(-.2,1,360),mean(Yall(chans,:,ltv.E<median(ltv.E)),3)')
subplot(212)
plot(linspace(-.2,1,360),mean(Yall(chans,:,ltv.E>median(ltv.E)),3)')

