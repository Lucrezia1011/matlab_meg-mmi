function [data,BadSamplesAll] = preproc_bids(data_name,highpass,lowpass,icaopt,plotopt)
% [data,BadSamples] = preproc_bids(data_name,highpass,lowpass,icaopt)
% icaopt = 1 to use ICA to remove eyeblinks and heartbeat
addpath('~/fieldtrip-20190812/fieldtrip_private')

% Data header info
hdr = ft_read_header(data_name);
% Get Bad channel names
fid = fopen([data_name,'/BadChannels']);
BadChannels = textscan(fid,'%s');
fclose(fid);

% get MEG channel names
channels = hdr.label(strcmp(hdr.chantype,'meggrad'));
% Delete Bad channels
chanInd = zeros(size(channels));
for iiC = 1:length(BadChannels{1})
    chanInd = chanInd | strcmp(channels,BadChannels{1}{iiC});
end
channels(find(chanInd)) = [];

%% Find large muscle artifacts
%  
% % With notch filter 
% cfg = [];
% cfg.dataset = data_name;
% cfg.continuous = 'yes';
% cfg.channel = channels;
% cfg.demean = 'yes';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [50 150];
% cfg.bsfilter = 'yes';
% cfg.bsfreq = [58 62; 118 122; 178 182];
% cfg.hilbert  = 'abs';
% data_muscle = ft_preprocessing(cfg);
% f = data_muscle.fsample;
% % find muscle artefact > 4pT (arbitrary threshold)
% % Include 0.5s either side of the artefact.
% data_muscle.trial{1} = data_muscle.trial{1} - mean(data_muscle.trial{1},2);
% [sampleMsl,~] = find(data_muscle.trial{1}'>4e-12);
% sampleMsl = unique(sampleMsl);
% warning('Found %d muscle artefacts/n',length(sampleMsl))
% BadMsl = zeros(length(sampleMsl),f+1);
% for iiM = 1:length(sampleMsl)
%     BadMsl(iiM,:) = sampleMsl(iiM)+[-0.5*f:0.5*f];
% end
% 
% BadMsl = unique(BadMsl(:));
% BadMsl(BadMsl<1) = []; % no negative samples

%%
 
% With notch filter 
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = channels;
cfg.demean = 'yes';
cfg.detrend = 'no';
cfg.bpfilter = 'no';
cfg.bpfreq = [0.5 300];
cfg.bsfilter = 'yes';
cfg.bsfreq = [58 62; 118 122; 178 182];

data = ft_preprocessing(cfg);
f = data.fsample;
% Eliminate bad trials(extract from ft_read_event)
[condNumbers,condLabels] = read_ctf_cls([data_name,'/ClassFile.cls']);

if any(strcmp(condLabels,'BAD'))
    condBad = condNumbers{strcmp(condLabels,'BAD')};
else
    condBad = [];
end
if any(strcmp(condLabels,'BAD_HeadMotion'))
    condBadHead =  condNumbers{strcmp(condLabels,'BAD_HeadMotion')};
else
    condBadHead = [];
end
% Find Trials marked as BAD
BadNumbers = unique([condBad, condBadHead]);

if ~isempty(BadNumbers) 
    BadTrials = zeros(data.hdr.nSamples,length(BadNumbers));
    for jN=1:length(BadNumbers)        
        BadTrials(:,jN) = (BadNumbers(jN)-1)*data.hdr.nSamples + (1:data.hdr.nSamples) ;       
    end
end
BadTrials = BadTrials(:);

% Find Bad data segments 
bsegFile =  [data_name,'/bad.segments'];
fid = fopen(bsegFile);
bsegs = cell2mat(textscan(fid,'%f'));
fclose(fid);

bsegs = reshape(bsegs,3,size(bsegs,1)/3)';
BadSegs = cell(size(bsegs,1),1);
if ~isempty(bsegs)   
    for jN=1:length(BadSegs)        
        BadSegs{jN} = (bsegs(jN,1)-1)*data.hdr.nSamples +...
            (round(bsegs(jN,2)*f):round(bsegs(jN,3)*f)) +1 ;       
    end  
end
BadSegs = cell2mat(BadSegs');
% Combine Bad Trials and Bad Segments
BadSamplesAll = unique([BadTrials; BadSegs']);
% Find Bad tail of dataset
indLast = find(diff(BadSamplesAll)~=1);
if isempty(indLast)
    indLast = 0;
end
BadSamplesLast = BadSamplesAll(indLast(end)+1:end);
BadSamples = BadSamplesAll(1:indLast(end));

% Eliminate Bad Segments at end of dataset
data.time{1}(BadSamplesLast) = [];
data.trial{1}(:,BadSamplesLast) = [];
data.sampleinfo = [1 length(data.time{1})];

% df = diff(data.trial{1},1,2);
% figure; histogram(df(:))
% % distribution of consecutive time points has stdev = 73.17 fT
%%

% % Baseline correct data by finding discontinuities from SQUID jumps 
% [sampleJump,sensJump] = find(abs(diff(data.trial{1}'))>5e-12); %jumps exceeding 1pT
% iia = ismember(sampleJump,BadSamples); % check if jumps are marked as Bad
% sampleJump(iia,:) = [];
% sensJump(iia,:) = [];

% Find large SQUID jumps
[sampleJump,sensJump] = find(abs(diff(data.trial{1}'))>100e-12); 

BadSQUID = cell(0);
iiN = 0;
for iib = 1:length(indLast)-1
    if iib == 1
        BadS = BadSamples(1):BadSamples(indLast(1));
    else
        BadS = BadSamples(indLast(iib-1)+1):BadSamples(indLast(iib));
    end
    iia = ismember(sampleJump,BadS);
    if nnz(iia) > 0
        iiN = iiN+1;
        BadSQUID{iiN} = BadS;
    end
    
end


% Check jumps after linear detrend: necessary for small jumps only
% indN = false(size(sampleJump));
% for iS = 1:length(sensJump)
%     tt = 1;
%     if sampleJump(iS)<=f*tt
%         s1 = 1:(sampleJump(iS)-1);
%         s2 = sampleJump(iS)+(1:f*tt);
%     elseif (sampleJump(iS)+f*tt)>data.sampleinfo(2)
%         s1 = sampleJump(iS)+(-f*tt:-1);
%         s2 = (sampleJump(iS)+1):data.sampleinfo(2);
%     else
%         s1 = sampleJump(iS)+(-f*tt:-1);
%         s2 = sampleJump(iS)+(1:f*tt);
%     end
%     % detrend
%     p = polyfit([s1,s2],data.trial{1}(sensJump(iS),[s1,s2]),1);
%     y1 = data.trial{1}(sensJump(iS),s1) - (p(1)*s1 + p(2));
%     y2 = data.trial{1}(sensJump(iS),s2) - (p(1)*s2 + p(2));
%    
%     indN(iS) = abs(mean(y1)-mean(y2))>2e-12;   
%     if indN(iS) && plotopt
%         figure
%         plot(s1(1):s2(end),data.trial{1}(sensJump(iS),s1(1):s2(end)))
%         hold on
%         plot([s1,s2],[y1,y2])
%         hold on
%         plot(sampleJump(iS),[-1 1]*5e-13)
%     end
% end
% 
% sampleJump = sampleJump(indN);
% sensJump = sensJump(indN);
% sensJumpu = unique(sensJump);
% 
% for iS = 1:length(sensJumpu)
%     sampleJumpi = sampleJump(sensJump==sensJumpu(iS));
%     
%     for iJ = 1:length(sampleJumpi)+1
%         if iJ==1
%             % First segment
%             s1 = 1;
%             s2 = sampleJumpi(iJ);
%             % mean correct
%             m1 = mean(data.trial{1}(sensJumpu(iS),s1:s2));
%             data.trial{1}(sensJumpu(iS),s1:s2) = ...
%                 data.trial{1}(sensJumpu(iS),s1:s2)-m1;
%             
%         else
%             % Start next segment where previous one ended
%             s1 = sampleJumpi(iJ-1)+1;   
%             if iJ > length(sampleJumpi)
%                 s2 = data.sampleinfo(2);
%             else
%                 s2 = sampleJumpi(iJ);
%             end
%             m1 = diff(data.trial{1}(sensJumpu(iS),[s1-1,s1]));
%             data.trial{1}(sensJumpu(iS),s1:s2) = ...
%                 data.trial{1}(sensJumpu(iS),s1:s2)-m1;
%         end
%         
%     end
% end

sensJumpu = unique(sensJump);
% separate channels with SQUID jumps
dataBadSQUID = data.trial{1}(sensJumpu,:);
% delete jump data segments
dataBadSQUID(:,cell2mat(BadSQUID)) = [];

% Filter data
data.trial{1} = data.trial{1} - mean(data.trial{1},2);
% Are we introducing filtering artefacts by applying highpass filter on
% data with discontinuities (i.e. after deleting bad segments)? Yes
% Apply highpass filter after correcting jumps, but before introducing
% discontinuities from bad data segments
data.trial{1} = ft_preproc_bandpassfilter(data.trial{1}, f, [highpass lowpass], [], [], [], []);

if ~isempty(sensJumpu)
    % Find clean jump
    [sampleJumpi,~] = find(abs(diff(dataBadSQUID'))>100e-12);
    if ~isempty(sampleJumpi)
    % correct jump
    for iS = 1:size(dataBadSQUID,1)
        
        for iJ = 1:length(sampleJumpi)+1
            if iJ==1
                % First segment
                s1 = 1;
                s2 = sampleJumpi(iJ);
                % mean correct
                m1 = mean(dataBadSQUID((iS),s1:s2));
                dataBadSQUID((iS),s1:s2) = ...
                    dataBadSQUID((iS),s1:s2)-m1;
            else
                % Start next segment where previous one ended
                s1 = sampleJumpi(iJ-1)+1;
                if iJ > length(sampleJumpi)
                    s2 = size(dataBadSQUID,2);
                else
                    s2 = sampleJumpi(iJ);
                end
                m1 = diff(dataBadSQUID((iS),[s1-1,s1]));
                dataBadSQUID((iS),s1:s2) = ...
                    dataBadSQUID((iS),s1:s2)-m1;
            end
            
        end
    end
    end
    % filter data
    dataBadSQUID = dataBadSQUID - mean(dataBadSQUID,2);
    dataBadSQUID = ft_preproc_bandpassfilter(dataBadSQUID, f, [highpass lowpass], [], [], [], []);
    % pad data channels with deleted jumps and substite into data struct
    indN = true(1,data.sampleinfo(2));
    for iS = 1:length(BadSQUID)
        indN(BadSQUID{iS}) = false;
    end
    dataBadSQUIDi = zeros(length(sensJumpu),data.sampleinfo(2));
    dataBadSQUIDi(:,indN) = dataBadSQUID;
    
    warning('Found %d SQUID jumps.\n',iiN)
    
    data.trial{1}(sensJumpu,:) = dataBadSQUIDi;
    clear dataBadSQUID dataBadSQUIDi
end

% Eliminate Bad Segments
data.time{1}(BadSamples) = [];
data.trial{1}(:,BadSamples) = [];
data.sampleinfo = [1 length(data.time{1})];

if plotopt == 1
    figure
    plot(data.trial{1}')
end

%%
% if ~isempty(sampleJump) 
%     keyboard %dbcont
% end
if icaopt == 1
    subN = data_name(1:9);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/',subN,'/',data_name(1:end-3)];
    if ~exist(data_path,'dir')
        mkdir(data_path)
    end
    
    rmpath(('~/fieldtrip-20190812/fieldtrip_private'))
    if exist([data_path,'/ICA_artifacts.mat'],'file')
        load([data_path,'/ICA_artifacts.mat']);
    else
        cfg = [];
        cfg.dataset = data_name;
        cfg.continuous = 'yes';
        cfg.channel = 'EEG';
        % cfg.demean = 'yes';
        % cfg.bpfilter = 'yes';
        % cfg.bpfreq = [1 150];
        try
            eog = ft_preprocessing(cfg);
            eog = eog.trial{1}(1,:);
        catch
            disp('Could not find EEG channel')
        end
        
        cfg =[];
        cfg.method = 'pca';
        comp_pca = ft_componentanalysis(cfg, data);
        score = comp_pca.trial{1}';
        compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
        icomp = nnz(compvar95) ;
        clc
        fprintf('%d components for 95perc. of data variance\n',icomp)
        
        if icomp>30
            disp('Reducing ICA components to 30')
            icomp = 30;
        end
        cfg =[];
        cfg.method = 'fastica';
        cfg.fastica.numOfIC = icomp;
        cfg.fastica.maxNumIterations = 100;
        comp = ft_componentanalysis(cfg, data);
        icomp = length(comp.label);
        
        figure
        cfg           = [];
        cfg.component = [1:icomp];       % specify the component(s) that should be plotted
        cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)
        
        
        cfg          = [];
        cfg.channel  = [1:5]; % components to be plotted
        cfg.viewmode = 'component';
        cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
        ft_databrowser(cfg, comp)
        
        if exist('eog','var')
            figure; % delete bad data segments
            eog(BadSamplesAll) = [];
            plot(abs(corr(eog',comp.trial{1}')),'*')
            grid on
            xlabel('ICA component')
            ylabel('Correlation with EOG')
        end
        icadel = input('ICA component to eliminate (input as [''01'';''02'']): ');
        
        cfg = [];
        cfg.channel = cell(size(icadel,1),1);
        for ii = 1:size(icadel,1)
            cfg.channel{ii}  = ['fastica0',icadel(ii,:)];
        end
        
        [comps] = ft_selectdata(cfg, comp);
        save([data_path,'/ICA_artifacts'],'comps')
        
    end
    close all
    cfg           = [];
    cfg.component = 1:length(comps.label);
    data          = ft_rejectcomponent(cfg, comps,data);
end
