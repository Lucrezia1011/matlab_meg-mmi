addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

meginfo = readtable('~/MEG_participantsinfo.csv');
%%
clc
iiS = 48;
sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan,'/meg'])

% bids naming convention
% MEG data:  sub-<label>_task-<label>_run-<index>_meg.ds
% MRI data:  sub-<label>_T1w
% behavioral :  sub-<label>_task-<label>_run-<index>_events.tsv

mrNames = dir;
mrNames = {mrNames.name};
ind = find(contains(mrNames,'.ds'));
%%
nInd = 3;
cfg=[];

cfg.method                  = 'copy';
% cfg.proc                    = 'raw';
cfg.dataset                 = mrNames{ind(nInd)};
cfg.bidsroot                = '/data/MBDU/bids/meg_mmi3/'; %top level directory for the BIDS output
cfg.sub                     = sdan; %subject name
cfg.run                     = 1; %optional
cfg.datatype                = 'meg'; %can be any of 'FLAIR', 'FLASH', 'PD', 'PDT2', 'PDmap', 'T1map', 'T1rho', 'T1w', 'T2map', 'T2star', 'T2w', 'angio', 'bold', 'bval', 'bvec', 'channels', 'coordsystem', 'defacemask', 'dwi', 'eeg', 'epi', 'events', 'fieldmap', 'headshape', 'ieeg', 'inplaneT1', 'inplaneT2', 'magnitude', 'magnitude1', 'magnitude2', 'meg', 'phase1', 'phase2', 'phasediff', 'photo', 'physio', 'sbref', 'stim'
cfg.participants.age        = meginfo.Age(iiS);
cfg.participants.sex        = meginfo.Sex(iiS); %'m' or 'f'
cfg.participants.group      = meginfo.Group(iiS);
% cfg.outputfile              = cfg.dataset;
% if nInd == 1
%     cfg.task                    = 'mmi3'; %task name is required for functional data
%     cfg.TaskName                = cfg.task;
%     cfg.TaskDescription         = 'Mood Machine Interface, 3-levels: high-low-high mood targets';
% else
%     cfg.task                    = 'rest';%'mmi3'; %task name is required for functional data
%     cfg.TaskName                = cfg.task;
%     cfg.TaskDescription         = 'Eyes open resting state scan.';
% end

if nInd < 3
    cfg.task                    = 'mmi3'; %task name is required for functional data
    cfg.TaskName                = cfg.task;
    cfg.TaskDescription         = 'Mood Machine Interface, 3-levels: high-low-high mood targets';
else
    cfg.task                    = 'rest';%'mmi3'; %task name is required for functional data
    cfg.TaskName                = cfg.task;
    cfg.TaskDescription         = 'Eyes open resting state scan.';
end

cd(cfg.dataset)
fileID = fopen([cfg.dataset(1:end-2),'infods'],'r');
TaskDate = [];
while isempty(TaskDate)
    tline = fscanf(fileID,'%s',1);
%     tline = fgetl(fileID);
    if contains(tline,'DATASET_COLLECTIONDATETIME')
        tline = fscanf(fileID,'%s',1);
        
        ind20 = strfind(tline,'20'); % find start of date, i.e. 2019 or 2020
        TaskDate = tline(ind20(1)+[0:13]);
    end
end
fclose(fileID);


cfg.scans.acq_time          =  [TaskDate(1:4),'-',TaskDate(5:6),'-',TaskDate(7:8),...
    'T',TaskDate(9:10),':',TaskDate(11:12),':',TaskDate(13:14)];
cd ..
disp(cfg.scans.acq_time)

data2bids(cfg);


% eval(sprintf(['!newDs -includeBad -includeBadChannels -includeBadSegments ',...
%     '%s /data/MBDU/bids/meg_mmi3/sub-%s/meg/%s'],mrNames{ind(1)},sdan))

%%
for iiS = 48
sdan = num2str(meginfo.SDAN(iiS));

cd(['/data/MBDU/bids/meg_mmi3/sub-',sdan])
if ~exist('anat','dir')
    mkdir anat 
    mkdir meg
    !mv *.ds meg
    !mv *.nii *.json anat
end

cd anat
% convert anatomical to bids
mrNames = dir;
mrNames = {mrNames.name};
ind = contains(mrNames,'.nii');
if nnz(ind)>1
    keyboard %dbcont
end
% Easier to rename than using datas2bids
% cfg =[];
% cfg.dataset                 = mrNames{ind};
% cfg.sub                     = sdan;
% cfg.method                  = 'copy';
% cfg.bidsroot                = '/data/MBDU/bids/meg_mmi3/';
% cfg.datatype                = 'T1w';
% cfg.acq                     = 'mprage';
% cfg.writejson               = 'no';
% data2bids(cfg);

indEx = strfind(mrNames{ind},'.');
mrEx = mrNames{ind}(indEx(1):end);  % file extension (.nii or .nii.gz)
eval(sprintf('!mv %s sub-%s_acq-mprage_T1w%s',mrNames{ind},sdan,mrEx))

ind = contains(mrNames,'.json');
eval(sprintf('!mv %s sub-%s_acq-mprage_T1w.json',mrNames{ind},sdan))
end