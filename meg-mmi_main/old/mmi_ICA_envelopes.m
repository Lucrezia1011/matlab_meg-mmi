clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

% Can do the same for ERF!
for sn = [1:9,11,12,14:16]% co-register 2 then run these! %[1,3,4,6,7,8,9,11,14,15,16]   %Subjects showing enough variation in mood
    
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
    clearvars -excephist subn sn data_names sub data_path runs freq_band freq
gridres = 8; % resolution of beamformer grid in mm
mniopt = true; % Use warped MNI grid to allow combination over sujects 
icaopt = false; % Use ICA to clean artefacts
dataprep = mmi_dataprep_mni(data_names{runs},gridres,mniopt,icaopt);
%%
data = dataprep.data;

cfg = [];
cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
data = ft_resampledata(cfg, data);

freq_band = [4 8; 8 13; 13 30; 1 35; 35 90];
for freq = [1:3,5]
filt_order = []; % default
data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, freq_band(freq,:),filt_order,'but');
% data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 ); 

data.trial{1} = data_filt;
clear data_filt

%% Read events
% 
% bv_match = match_triggers_fc(dataprep.dataname, 300);
% 
% answer_match = bv_match.answer;
% choice_match =  bv_match.choice;
% outcome_match  = bv_match.outcome;
% mood_match = bv_match.ratemood;
% blockmood_match = bv_match.blockmood;
% slider_match = bv_match.slider;
% blockslider_match = bv_match.blockslider;
% ITI_match = bv_match.ITI ;
% buttonpress = bv_match.buttonpress;
% 
% pRPE = outcome_match.win == 1 ;
% nRPE = outcome_match.win == -1 ;
% pRPE_sample = outcome_match.sample(pRPE);
% nRPE_sample = outcome_match.sample(nRPE);


%% Beamfomer
icacomps = length(dataprep.data.cfg.component);

C = cov(data.trial{1}');
E = svd(C); 
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

Cr = C + 4*noiseC; % need to normalise because of ICA

L = dataprep.leadfield.leadfield;
VE(1:size(L,2)) = {0};
W(1:size(L,2)) = {0};
%
% 
% sdata = svd(data.trial{1}');
% noise = sdata(end-icacomps);

parfor ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)       
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
%         if d(3) < 1
%             jj = 2; % The minumum singular value is degenerate
%         else
%             jj =3;
%         end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        w = Cr\lfo / (lfo'/Cr*lfo) ;             
%         w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;       
        W{ii} = w;
    end
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end    
end
clc
fprintf('Beamformer finished\n' )

for ii = 1:length(L)
    if ~isempty(L{ii})
        w = W{ii};
%         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
%         Better Hall's or normalized weights
        wnorm = w/sqrt(w'*noiseC*w);    
        VE{ii} =  wnorm'*data.trial{1};
    end
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
end

%% Calculate hilbert envelope and downsample to 1Hz

nsamples = length(data.time{1})/data.fsample; % Number of samples at 1Hz
nvox  = nnz(dataprep.leadfield.inside);

VEhilb(1:size(L,2)) = {0};
% VEhilb = zeros(nsamples, nvox);
f = data.fsample;
for ii = 1:length(L)
    if ~isempty(L{ii})
        
        h = abs(hilbert(VE{ii}));
%         VEhilb{ii} = downsample(h,f);
        VEhilb{ii} = resample(h,1,f);
    end
    if mod(ii,300) == 0
        fprintf('Calculate dowsampled hilbert envelopes %.1f perc.\n', ii/length(L)*100 )
    end
end

%%
% sourcemodel = dataprep.sourcemodel;
save([data_path,data_names{runs},'/beamforming/VE_ICA_envelopes_',...
    num2str(freq_band(freq,1)),'-',num2str(freq_band(freq,2)),'Hz.mat'],'VEhilb')
end
end
end
%% Load MNI brain and source array
clearvars -except subn gridres

ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));

spmpath = '/spin1/home/linux/liuzzil2/fieldtrip-20190812/external/spm8';

% mri = ft_read_mri([spmpath,'/templates/T1.nii']);
mri = ft_read_mri('~/MNI152_T1_2009c.nii');

nlocs = length(sourcemodel.inside);

%% Combine time courses from all subjects

% Number of datasets per subject
nsets = [1,1,1,1,1,1,3,1,1,1,2,1,1,1,2,2];
% subjects analyzed
slist = [1:9,11,12,14:16]; 

VE = cell(sum(nsets(slist)),nlocs);
ii = 0;

for sn = slist %[1,3,4,6,7,8,9,11,14,15,16] %already ran n.6   %Subjects showing enough variation in mood
    
    sub = subn(sn,:);
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    cd(data_path)
    name_list = dir;
    data_names = cell(1);
    for k = 1:length(name_list)
        if strncmp(name_list(k).name, data_name, 18) 
            ii = ii+1;
            load([data_path,name_list(k).name,'/beamforming/VE_ICA_envelopes_13-30Hz.mat'])
%             delete([data_path,name_list(k).name,'/beamforming/VE_ICA_enveopes.mat'])
            VE(ii,:) = VEhilb;
            clear VEhilb
        end
    end

end


