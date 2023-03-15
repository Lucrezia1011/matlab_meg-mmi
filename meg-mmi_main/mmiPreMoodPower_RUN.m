% Lucrezia Liuzzi, last updated 2021/03/18
%
% Run mmiPreMoodPower (if needed) and combine results from all subjects to 
% save for linear mixed model analysis
% select exploratory, confirmatory or post-hoc analysis
%
% mmiPreMoodPower(data_name,roiopt,gridres,freqband,mu)
% data_name = name of dataset (.ds)
% roiopt    = 'grid' beamformer on mni grid, 'sens' sensor level
% gridres   = grid resolution in mm (for beamformer)
% freqband  = frequency band [low_f, high_f]
% mu        = beamformer regularization parameter, e.g. mu=0.05 (fraction of maximum singular value of covariance)


clear all
close all
% clc
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

%% Run mmiPreMoodPower
addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


gridres= 5;
% mu = 0.002; %saved as mu 0 
mu = 0.05; %mu=0.01

freqband = [25 40]; 

for ii = 1:length(data_list)
    data_name = data_list{ii};
    sub = data_name(5:9);
    mmiPreMoodPower(data_name,roiopt,gridres,freqband,mu); % mu = 0.05 
%     mmiPreMoodPower_multiSpheres(data_name,roiopt,gridres,freqband,mu); % mu = 0.05 
end

return

%% Find all common gradiometers for sensor space analysis

BadChannels = [];
for s = 1:length(data_list)
    
    data_name = data_list{s};
    sub = data_name(5:9);
    cd(['/data/MBDU/MEG_MMI3/data/bids/sub-',sub,'/meg/'])
    
    % Get Bad channel names
    fid = fopen([data_name,'/BadChannels']);
    BadChannel = textscan(fid,'%s');
    fclose(fid);
    BadChannels = cat(1,BadChannels,BadChannel{1});
    
end

BadChannels = unique(BadChannels);
% get MEG channel names
hdr = ft_read_header(data_name);
channels = hdr.label(strcmp(hdr.chantype,'meggrad'));
% Delete Bad channels
chanInd = zeros(size(channels));
for iiC = 1:length(BadChannels)
    chanInd = chanInd | strcmp(channels,BadChannels{iiC});
end
channels(chanInd) = [];

% save('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/sensors','channels'); 
save('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/confirm/sensors','channels');
%% Combine data from all subjects
Pall = [];
Vall = [];

ltvMood = [];

Ytfs = [];
Ytfsp = [];

if strcmp(roiopt,'sens')

    if ~exist('sensall','var')
        load('/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/confirm/sensors.mat')
    end

    cd(['/data/MBDU/MEG_MMI3/data/bids/',data_list{1}(1:9),'/meg/'])
    hdr = ft_read_header(data_list{1});
    channelsall = hdr.label(strcmp(hdr.chantype,'meggrad'));
end

for sn = 1:length(data_list) 
    
    
    data_name = data_list{sn};
    sub = data_name(5:9);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/',data_name(1:end-3)];
    
    
    if strcmp(roiopt,'grid')
        cd(data_path)
        load('leadfields_5mm.mat')
        if exist('gridall','var')
            gridall = gridall & grid.inside;
        else
            gridall = grid.inside;
        end
        
%         load(sprintf('pre_mood_grid_%.0f-%.0fHz_mu%g',...
%             freqband(1),freqband(2),mu*100))
        load(sprintf('pre_mood_grid_%.0f-%.0fHz',...
            freqband(1),freqband(2)))
        Pmni = zeros(size(grid.inside,1),size(P,2));
        Pmni(grid.inside,:) = P;
        Pall = cat(2,Pall,Pmni);
        
        
    elseif  strcmp(roiopt,'sens')
        
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
        cd(data_path)
        
        load(sprintf('pre_mood_sens_%.0f-%.0fHz',...
            freqband(1),freqband(2)))
        
        [~,~,ind]= intersect(channels,channelSub);
        Vall = cat(2,Vall,V(ind,:));
   
    end
    
    ltvmood.recording = repmat(sn,size(ltvmood,1),1);
    if isempty(ltvMood)
        ltvMood = ltvmood;
    else
        ltvMood(end+(1:size(ltvmood,1)),:) = ltvmood;
    end
    
    
end

clear ltvmood TFS P Pmni grid
if strcmp(roiopt,'grid')
    Pall = Pall(gridall,:); 
end

%% Write data

if strcmp(roiopt,'grid')
    % Write data (grid)
    out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_grid/pre_mood/confirm/';
    % dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz.txt',out_path,...
    %     freqband(1),freqband(2)),Pall);
    dlmwrite(sprintf('%s/powergrid_%.0f-%.0fHz_mu%g.txt',out_path,...
        freqband(1),freqband(2),mu*100),Pall);
    dlmwrite([out_path,'/mni_grid.txt'],gridall);
    writetable(ltvMood,[out_path,'/latent_vars.csv']);

elseif strcmp(roiopt,'sens')
    % Write data (sens)
    out_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sens/pre_mood/confirm/';
    dlmwrite(sprintf('%s/powersens_%.0f-%.0fHz.txt',out_path,...
        freqband(1),freqband(2)),Vall);
    writetable(ltvMood,[out_path,'/latent_vars.csv']);
end
