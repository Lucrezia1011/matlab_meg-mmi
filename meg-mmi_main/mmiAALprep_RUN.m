% Lucrezia Liuzzi, last updated 2021/03/15
%
% Calculate evoked responses (lowpass 30Hz) to gamble feedback
% Saves timecourses and corresponding mood model parameters as .mat file
%
% mmiAALprep(data_name,twind)
% data_name = name of dataset (.ds)
% twind     = time window in seconds [t1, t2]
% roiopt    = 'AAL' or 'sens', beamform on AAL atlas or keep sensor array
% opt       = cell array with trigger selection, e.g. {'outcome';'cue'}
%             'outcome': gamble feedback
%             'cue'    : gamble options presentation
%             'choice' : gamble choice selection
% Warning: data path and output directory are hard-coded!

clear all
close all
clc
meginfo = readtable('~/MEG_participantsinfo.csv'); % on github
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
subexclude = [10,24]; % sub-24199, sub-24128


analy_case = 'confirm';
roiopt = 'grid'; % running for grid
switch roiopt
    case 'grid'
        subexclude = [subexclude,26,49,53];
end

switch analy_case
    case 'explore'
        subexclude = [subexclude,13,17:56];  % exploratory 14 subjects
    case 'confirm'
        subexclude = [subexclude,1:12,14:16]; % confirmatory 37 subjects
    case 'posthoc'
       % all available subjects  
end

roiopt = 'AAL'; 
switch roiopt
    case 'AAL'
        subexclude = [subexclude,26,49,53]; % 3 subjects with no mri co-registration
end

Nlist(subexclude) = []; 
zz= 0;
for sn = Nlist 
        
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

addpath /home/liuzzil2/fieldtrip-20190812/ % fieldtrip path
ft_defaults

twind = [-.2,1];
condition = 'outcome';% 'outcome', 'cue' , 'choice'

% Run mmiAALprep 
% for ii = 1:length(data_list) 
%     data_name = data_list{ii};
%     sub = data_name(5:9);
%     processing_folder = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3),'/'];
%     mmiAALprep(data_list{ii},twind,roiopt,condition)
% 
% end
% return
%%

derivative_dir = '/data/MBDU/MEG_MMI3/data/derivatives/';
r = 0;
Yall = [];
ltv = [];

Ym = [];

ntrialsOut = [];

if strcmp(roiopt,'sens')
    load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors.mat')
    cd('/data/MBDU/MEG_MMI3/data/bids/sub-24071/meg/')
    hdr = ft_read_header('sub-24071_task-mmi3_run-1_meg.ds');
    channelsall = hdr.label(strcmp(hdr.chantype,'meggrad'));
end
time = linspace(twind(1),twind(2),360);

for sn = 1:length(data_list) 
    
    data_name = data_list{sn};   
    sub = data_name(5:9);
    processing_folder = [derivative_dir,'sub-',sub,'/',data_name(1:end-3),'/'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save_name = [processing_folder,'evoked_',condition,'_',roiopt,'.mat'];
    load(save_name)
    
    switch condition
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
        
     end
    
    
    Yall = cat(3,Yall,Yout);
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
clear ltvchoice ltvcue ltvout Y Ytheta Ycue Yout Ychoice
%% Sign flip for source localized evoked responses
if strcmp(roiopt,'AAL')

    Sout = zeros(size(Yall));

    for ii = 1:size(Ym,1)
        Yu = squeeze(Ym(ii,:,:))';
        Yu = zscore(Yu,0,2);

        % first pass
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

    figure; set(gcf,'color','w')
    subplot(131); imagesc(mean(Yall,3)); caxis([-1 1]*.3)
    title(sprintf('Grandaverage over all subjects\n of beamformed evoked responses'))
    xlabel('timepoint'); ylabel('ROI')

    Yall = Yall.*Sout;

    subplot(132); imagesc(mean(Yall,3)); caxis([-1 1]*.3)
    title(sprintf('Maximised ROI temporal\n correlation over subject'))
    xlabel('timepoint'); ylabel('ROI')
    % Flip ROIs to match each other
    % Yu = cat(2,mean(YCue,3),mean(YChoice,3),mean(YOut,3));
    Yu = mean(Yall,3);

    % second pass
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
    Yall = Yall.*Sout;

    subplot(133); imagesc(mean(Yall,3));  caxis([-1 1]*.3)
    title(sprintf('Maximised temporal correlation\n over all ROIs'))
    xlabel('timepoint'); ylabel('ROI'); 
    c = colorbar;
    c.Position = [0.93 c.Position(2:4)];
end

%% Write data

switch roiopt
    case 'AAL'
        if strcmp(condition,'outcome')
            
            Yall = permute(Yall,[2,1,3]);
            Yall = reshape(Yall,size(Yall,1)*size(Yall,2),size(Yall,3)); % nrois * npoints * ntrials

            dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/meg_trials_evoked_',condition,'.txt'],Yall);
            writetable(ltv,['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/confirm/latent_vars_evoked_',condition,'.csv']);
        else
            subs = unique(ltv.subject);
            Ym = zeros(size(Yall,1),size(Yall,2),length(subs));
            for s = 1:length(subs)
                Ym(:,:,s) = mean(Yall(:,:,ltv.subject == subs(s)),3);
            end
            Ym = permute(Ym,[2,1,3]);
            Ym = reshape(Ym,size(Ym,1)*size(Ym,2),size(Ym,3)); % nrois * npoints * ntrials

            dlmwrite(['/data/MBDU/MEG_MMI3/results/mmiTrial_aal/meg_trials_evoked_',condition,'.txt'],Ym);
        end
    case 'sens'
        writetable(ltv,['/data/MBDU/MEG_MMI3/results/mmiTrial_sens/evoked_outcome/latent_vars_evoked_',condition,'.csv']);
        
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
