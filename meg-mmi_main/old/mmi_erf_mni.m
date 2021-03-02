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
for sn = [14:16,7,3,12] %[7,11,12,15,16]  % co-register 2 then run these! %[1,3,4,6,7,8,9,11,14,15,16]   %Subjects showing enough variation in mood
    
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
    clearvars -except subn sn data_names sub data_path runs freq_band freq
gridres = 5; % resolution of beamformer grid in mm
mniopt = true; % Use warped MNI grid to allow combination over sujects 
icaopt = true;
dataprep = mmi_dataprep_mni(data_names{runs},gridres,mniopt,icaopt);
%%
data = dataprep.data;

% cfg = [];
% cfg.resamplefs = 300; % Downsample to 300Hz for ease of memory
% data = ft_resampledata(cfg, data);


filt_order = []; % default
% data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, freq_band(freq,:),filt_order,'but');
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 ); 

data.trial{1} = data_filt;
clear data_filt

%% Read events

bv_match = match_triggers_fc(dataprep.dataname);

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;
buttonpress = bv_match.buttonpress;
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

% Cr = C + 4*noiseC; % old normalization
Cr = C + 0.01*eye(nchans)*E(1); % 1% max singular value

L = dataprep.leadfield.leadfield;

VEp(1:size(L,2)) = {0};
VEn(1:size(L,2)) = {0};
% VEg(1:size(L,2)) = {0};
% VEs(1:size(L,2)) = {0};
W(1:size(L,2)) = {0};

% dataerf = define_trials(buttonpress,data,bv_match.time,[-0.2 1]);

Rsamples = outcome_match.sample(outcome_match.win~=0);
% dataerf = define_trials(Rsamples,data,bv_match.time,[-0.2 1]);
[dataerf, ttdel] = define_trials(Rsamples,data,bv_match.time,[-1 2]); % new window lenght


sdata = svd(cell2mat(dataerf.trial)');
noise = sdata(end-icacomps);


ntrials = length(dataerf.time);
nsamples = length(dataerf.time{1});

dataerfc = cell2mat(dataerf.trial);

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
        
%         wnorm = w/sqrt( sum( (w*noise).^2) ); % Not good normalization!
%         Better Hall's or normalized weights
        wnorm = w/sqrt(w'*noiseC*w);    
%         w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;       
        W{ii} = wnorm;
    end
    if mod(ii,1000) == 0
        clc
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end    
end
clc
fprintf('Beamformer finished\n' )

trialsRPE = outcome_match.RPE(outcome_match.win~=0);
trialsRPE(ttdel) = [];
Rwin = trialsRPE>0; 
Rgamble = abs(trialsRPE)>5; % big gambles!
Rsafe = abs(trialsRPE)<3; % safe gambles

for ii = 1:length(L)
    if ~isempty(L{ii})
        w = W{ii};

        dataloc = w'*dataerfc;
        dataloc = reshape(dataloc, [nsamples, ntrials]);
        
        % Added on 08/01/2020
%         dataloc = zscore(dataloc);
        VEp{ii} = mean(dataloc(:,Rwin),2); % All wins
        VEn{ii} = mean(dataloc(:,~Rwin),2); % All losses
%         VEg{ii} = mean(dataloc(:,Rgamble),2); % 
%         VEs{ii} = mean(dataloc(:,Rsafe),2);
 
    end
    if mod(ii,1000) == 0
        clc
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
end

%%

%%
% sourcemodel = dataprep.sourcemodel;
save([data_path,data_names{runs},'/beamforming/VE_erf_mu1.mat'],'VEp','VEn','Rwin')
end
end
