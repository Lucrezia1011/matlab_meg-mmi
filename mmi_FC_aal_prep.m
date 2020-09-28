function mmi_FC_aal_prep(data_name,twind,inducedopt,roiopt,Popt,FCopt)
% roiopt = 'g' guassian weighting
% roiopt = 'c' centroid

% addpath /home/liuzzil2/fieldtrip-20190812/
% ft_defaults
% addpath('~/fieldtrip-20190812/fieldtrip_private')
% addpath ~/ppyll1/matlab/svdandpca

%% Co-register MRI from fiducial positions

% LTA model latent variables:
% EC: Expectation of certain value
% EG: Expectation during gabling
% Ediff: Drift rate
% LTA: Long term average with gamma:   1/t * sum_i=1 ^t(V(i)^gamma),   cumsum(LTA.OutcomeAmount^gamma)./(1:ntrials)'
% V_i^gamma = outcome of trial i
% new_p = subjective winning probability
% RPE = Reward prediction error
% LTA_sum  = sum(LTA)
% RPE_sum = sum(RPE)
% log_like
% mood_log_like

sub = data_name(1:5);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)




%% Co-register MRI

processing_folder = [data_path,data_name,'/beamforming'];
if ~exist(processing_folder,'dir')
    mkdir(processing_folder)
end

mri_name = [sub,'_anat+orig.BRIK'];

if ~exist(mri_name,'file')
    unix(['gunzip ',mri_name])
end

mri = ft_read_mri(mri_name,'dataformat','afni_brik');

tagset_shape = mri.hdr.TAGSET_NUM;
tagset_coord = mri.hdr.TAGSET_FLOATS;
tagset_coord = reshape(tagset_coord,fliplr(tagset_shape)); % nas, lpa, rpa

tagset_p = zeros(1,3);  % Ideal orientation {RL; PA; IS}
for ii =1:3
    if strcmp(mri.hdr.Orientation(ii,:),'AP') || strcmp(mri.hdr.Orientation(ii,:),'PA')
        tagset_p(ii) = 2;
    elseif strcmp(mri.hdr.Orientation(ii,:),'LR') || strcmp(mri.hdr.Orientation(ii,:),'RL')
        tagset_p(ii) = 1;
    elseif strcmp(mri.hdr.Orientation(ii,:),'SI') || strcmp(mri.hdr.Orientation(ii,:),'IS')
        tagset_p(ii) = 3;
    end
end

m = [   -1  0   0   mri.dim(1)
    0   -1  0   mri.dim(2)
    0   0   1   1
    0   0   0   1] ;


tagset_coord = tagset_coord(tagset_p,:)'; % fiducials have shuffled coordinates

mri.transform(1:3,4) = mri.hdr.ORIGIN; % change translation to origin

mri.transform = mri.transform/m;
fiducial_coord = (mri.transform \[tagset_coord,ones(3,1)]')';

cfg = [];
cfg.method = 'fiducial';
cfg.fiducial.nas    = fiducial_coord(1,1:3); %position of nasion
cfg.fiducial.lpa    = fiducial_coord(2,1:3); %position of LPA
cfg.fiducial.rpa    = fiducial_coord(3,1:3); %position of RPA
cfg.coordsys = 'ctf';
cfg.viewresult = 'no';

mri = ft_volumerealign(cfg,mri);

if ~exist([sub,'_coreg.nii'],'file')
    writebrik([sub,'_coreg'],mri);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment MRI
if ~exist([processing_folder,'/headmodel.mat'],'file')
    cfg = [];
    cfg.output  = 'brain';
    segmentmri = ft_volumesegment(cfg,mri);
    
    % Head model
    
    cfg = [];
    cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, segmentmri);
    
    save([processing_folder,'/headmodel.mat'],'vol')
else
    load([processing_folder,'/headmodel.mat']);
end
sens = ft_read_sens(data_name,'senstype','meg');


%% AAL atlas
gridres = 5; % resolution of beamformer grid in mm

% Load fieldtrip 10mm MNI grid
ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));
template_grid = sourcemodel;
atlas = ft_read_atlas('~/fieldtrip-20190812/template/atlas/aal/ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,sourcemodel.unit);

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodelAAL = ft_sourceinterpolate(cfg, atlas, sourcemodel);

clear sourcemodel

%% Sourcemodel warp MNI grid

% sourcemodel based on 5mm grid MNI brain
cfg = [];
cfg.mri = mri;
cfg.warpmni = 'yes';
cfg.template  = template_grid; % Has to be template grid! Made from ft_prepare_sourcemodel
cfg.unit      = 'm';
cfg.nonlinear = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);

%% Find location of AAL ROIs
R = length(sourcemodelAAL.tissuelabel);
locs = zeros(R,3);
locsAAL = cell(R,1);
for ii = 1:R
    ind = find(sourcemodelAAL.tissue == ii);
    voxc = mean(sourcemodel.pos(ind,:)); % centroid
    locs(ii,:) = voxc;
    
    locsAAL{ii} = sourcemodel.pos(ind,:);
    
end

if strcmp(roiopt,'g')
    locsc = locs;
    locs = cell2mat(locsAAL);
end
%% Calculate lead fields

cfg                 = [];
cfg.grad            = sens;
cfg.headmodel       = vol;
cfg.reducerank      = 2;
cfg.channel         = {'MEG'};
cfg.sourcemodel.pos = locs;
cfg.sourcemodel.unit   = 'm';
cfg.siunits         = true;
cfg.normalize = 'no'; % To normalize power estimate (center of the head bias for beamformer and superficial bias for mne)
[grid] = ft_prepare_leadfield(cfg);

%% Clean data with ICA

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'MEG';
% cfg.demean = 'yes';
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [1 150];
data = ft_preprocessing(cfg);
f = data.fsample;

if exist([processing_folder,'/ICA_artifacts.mat'],'file')
    load([processing_folder,'/ICA_artifacts.mat']);
    
end

cfg           = [];
cfg.component = 1:length(comps.label);
data          = ft_rejectcomponent(cfg, comps,data);

%%
filt_order = []; % default
if strncmp(inducedopt{1},'gamma',5)
    data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[40 150],filt_order,'but');
else
    data_filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample,[1 40],filt_order,'but');
end

data.trial{1} = data_filt;
clear data_filt

%% Beamfomer
icacomps = length(data.cfg.component);

C = cov(data.trial{1}');
E = svd(C);
nchans = length(data.label);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

% Cr = C + 4*noiseC; % old normalization
Cr = C + 0.05*eye(nchans)*E(1); % 5% max singular value

if strcmp(roiopt,'g')
    VE = cell(R,1);
    n =0;
    for r = 1:R
        clc
        fprintf('SAM running %d/%d .\n', r,R)
        
        L = grid.leadfield( n + (1:size(locsAAL{r},1)) );
        
        VEr = zeros(data.sampleinfo(2),size(locsAAL{r},1));
        
        voxc = locsc(r,:); % centroid
        GD = zeros(1,size(locsAAL{r},1));
        for ii = 1:length(L)
            
            d = sqrt(sum((grid.pos(n+ii,:)-voxc).^2,2)); % distance from centroid
            GD(ii) = exp(-(d.^2)/1e-4); % gaussian weigthing
            lf = L{ii}; % Unit 1Am
            if GD(ii) > 0.05 && ~isempty(lf) % include voxels with weighting > 5%
                % %  G O'Neill method, equivalent to ft
                [v,d] = svd(lf'/Cr*lf);
                d = diag(d);
                jj = 2;
                
                lfo = lf*v(:,jj); % Lead field with selected orientation
                
                w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
                
                VEr(:,ii)  = GD(ii)*w'*data.trial{1};
                
            end
        end
        
        sf = corr(VEr); % check sign
        [~,ind] = max(GD);
        sf= sign(sf(ind,:));
        sf(isnan(sf)) = 0;
        VEr = VEr.*sf;
        VE{r} = sum(VEr,2);
        n = n + size(locsAAL{r},1);
    end
    
else
    L = grid.leadfield;
    
    VE(1:size(L,2)) = {0};
    W(1:size(L,2)) = {0};
    for ii = 1:length(L)
        lf = L{ii}; % Unit 1Am
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
        
        lfo = lf*v(:,jj); % Lead field with selected orientation
        
        w = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        W{ii} = w;
        VE{ii}  = w'*data.trial{1};
        
        clc
        fprintf('SAM running %d/%d .\n', ii, R)
        
    end
end
%     fprintf('Beamformer finished\n' )

VE = cell2mat(VE');
% zscore VE, check dimensions
% VE = zscore(VE);
fprintf('Done.\n')

%% Read events

[bv_match,~] = match_triggers_fc(data_name);

% cue_match = bv_match.answer;
% choice_match = bv_match.choice;
% outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
% blockmood_match = bv_match.blockmood;
tasktime = bv_match.time;

% ntrials = nnz(~isnan(bv.outcomeAmount));
% inds = find(~isnan(bv.outcomeAmount));
% ind1 = inds(1)-1;
%
% hsind = find(~isnan(bv.happySlider_response));
% mood_match.mood(hsind) =  bv.happySlider_response(hsind);
% Fsample = griddedInterpolant(find(mood_match.sample),mood_match.sample(mood_match.sample~=0),'linear');
% mood_match.sample(hsind) = Fsample(hsind);
%
% bv = bv(inds,:);
% bv.trialNumber = (1:ntrials)'-1;

%%

save_name = ['/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/FC/',sub];

n = str2double(data_name(end-3));
if ~isnan(n) %check for number at end of filename
    save_name = [save_name,'_',data_name(end-3)];
else
    save_name = [save_name,'_1'];
end


datave = data;
datave.trial{1} = VE';
datave.label = sourcemodelAAL.tissuelabel';
dataf = datave;

clear VE data



clear ltvchoice ltvout ltvcue Ychoice Yout Ycue Yalpha Ybeta Ytheta

for freq = 1:length(inducedopt)
    % Induced oscillations
    switch inducedopt{freq}
        case 'delta'
            freql = [1 4];
            filtertype = 'firls';
        case 'theta'
            freql = [4 8];
            filtertype = 'but';
        case 'alpha'
            freql = [8 13];
            filtertype = 'but';
        case 'beta'
            freql = [13 30];
            filtertype = 'but';
        case 'gamma'
            freql = [40 100];
            filtertype = 'but';    
    end
    
    
    data_filt = ft_preproc_bandpassfilter(datave.trial{1}, datave.fsample, freql,[],filtertype);
    
    %         data_lk = econleakagecorr(data_filt','y');
    %         % Run leakage correction multivariate
    %
    %         env = abs(hilbert(data_lk));
    % %         envd = resample(env,10,data.fsample);
    % %         figure;
    % %         imagesc(corr(envd)); colorbar; caxis([0, 0.2])
    %
    %         dataf.trial{1}= env';
    % %         dataf.trial{1}= zscore(data_filt)';
    %
    %         twinde = twind/2*[-1,1];
    %         [dataout,ttdel]= define_trials(mood_match.sample(mood_match.sample~=0), dataf, tasktime, twinde,0);
    %
    %         nchans = length(dataout.label);
    %         ntrials = length(dataout.trial);
    %         S = str2double(sub)*ones(ntrials,1);
    %
    % %         AEC = cell(1,ntrials);
    %         AECs = zeros(nchans,ntrials);
    %         for k = 1:ntrials
    %             envd = resample(dataout.trial{k}',1,dataout.fsample);
    %             aec = corr(envd);
    %             AECs(:,k) = (sum(aec)-1)/(nchans-1);
    %         end
    
    %%
    
    if Popt == 1
        
        dataf.trial{1}= abs(hilbert(data_filt'))';
        twinde = twind/2*[-1,1];
        [dataout,ttdel]= define_trials(mood_match.sample(mood_match.sample~=0), dataf, tasktime, twinde,0);
        
        nchans = length(dataout.label);
        ntrials = length(dataout.trial);
        S = str2double(sub)*ones(ntrials,1);
        
        Pow = zeros(nchans,ntrials);
        
        for k = 1:ntrials
            datak = dataout.trial{k};
            Pow(:,k) = mean(datak,2); % Average power within window
            
        end
        %%
        
        ltvind = mood_match.bv_index(mood_match.sample~=0); % indeces start at 13
        ltvind(ttdel) = [];
        
        mood = mood_match.mood(ltvind);
        trials = mood_match.bv_index(ltvind)-12;
        
        ltvmood = table(S,trials',mood','VariableNames',...
            {'subject','trial','mood'});
        
        %%
        
        eval(sprintf('Pow_%s = Pow; ',inducedopt{freq}));
        
        if exist([save_name,'.mat'],'file')
            save(save_name,'ltvmood',['Pow_',inducedopt{freq}],'-append');
        else
            save(save_name,'ltvmood',['Pow_',inducedopt{freq}]);
        end
        clear Pow
        
    end
    
    
    
    
    if FCopt == 1
        
        29,30 % right insula
        31,32 % ACC
        69,70 % paracentral lobule
        45,46 % Cuneus
        
        oind = outcome_match.sample(outcome_match.win==-1);
        tind = oind+(-.5*f:f*1)';      
        iac = zeros(f*1.5+1,2);
        
        for n = 1:2
            if n == 1
                ii =30; jj =45;
            else
                ii =45; jj =30;
            end
        x = VEz(:,ii);  y = VEz(:,jj);
        y = leakage_reduction(y,x);     
        x = abs(hilbert(x));
        y = abs(hilbert(y));
        x = x(tind); y = y(tind);      
        iac(:,n) = mean(x.*y,2);
        end     
        hold on
        plot(linspace(-.5,1,f*1.5+1),mean(iac,2))
        
        dataf.trial{1}= data_filt;
        twinde = twind/2*[-1,1];
        [dataout,ttdel]= define_trials(mood_match.sample(mood_match.sample~=0), dataf, tasktime, twinde,0);
        
        nchans = length(dataout.label);
        ntrials = length(dataout.trial);
        S = str2double(sub)*ones(ntrials,1);
        
        fsamp = 10;
        %         AEC = cell(1,ntrials);
        
        AECs = zeros(nchans,ntrials);
        inds = find(triu(ones(nchans),1));
        AECk = zeros(length(inds),ntrials);
        for k = 1:ntrials
            datak = dataout.trial{k}';
            aec = zeros(nchans);
            for ii = 1:nchans
                x = datak(:,ii);
                envx = abs(hilbert(x));
                envx([1:20,end-19:end]) = [];
                envx = resample(envx,fsamp,f);
                for jj = [1:ii-1,ii+1:nchans]
                    y = datak(:,jj);
                    y = leakage_reduction(y,x);
                    envy = abs(hilbert(y));
                    envy([1:20,end-19:end]) = [];
                    envy = resample(envy,fsamp,f);
                    aec(ii,jj) = corr(envx,envy);
                end
                clc; fprintf(['Subject %s, freq. %s\nTime segment %.0f/%.0f\n'...
                    'Done ROI %.0f\n'],sub,inducedopt{freq},k,ntrials,ii)
            end
            aec = (aec+aec')/2;
            AECs(:,k) = sum(aec)/(nchans-1);
            AECk(:,k) = aec(inds);
        end
        %%
        
        ltvind = mood_match.bv_index(mood_match.sample~=0); % indeces start at 13
        ltvind(ttdel) = [];
        
        mood = mood_match.mood(ltvind);
        trials = mood_match.bv_index(ltvind)-12;
        
        ltvmood = table(S,trials',mood','VariableNames',...
            {'subject','trial','mood'});
        
        %%
        
        eval(sprintf('AEC_%s = AECs; ',inducedopt{freq}));
        eval(sprintf('Taec_%s = AECk; ',inducedopt{freq}));
        if exist([save_name,'.mat'],'file')
            save(save_name,'ltvmood',['AEC_',inducedopt{freq}],['Taec_',inducedopt{freq}],'-append');
        else
            save(save_name,'ltvmood',['AEC_',inducedopt{freq}],['Taec_',inducedopt{freq}]);
        end
        clear AECs AECk aec
        
    end
end




