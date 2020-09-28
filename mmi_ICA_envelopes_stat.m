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
gridres =8;
%% Load MNI brain and source array

ftpath   = '/home/liuzzil2/fieldtrip-20190812/';
load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d',num2str(gridres),'mm']));

spmpath = '/spin1/home/linux/liuzzil2/fieldtrip-20190812/external/spm8';

% mri = ft_read_mri([spmpath,'/templates/T1.nii']);
mri = ft_read_mri('~/MNI152_T1_2009c.nii');

nlocs = length(sourcemodel.inside);

%% Combine time courses from all subjects

% Number of datasets per subject
nsets = [1,1,1,1,1,1,2,1,1,1,2,1,1,1,2,2];
% subjects analyzed
slist = [1:9,11,14:16]; %12 is too short

VE = cell(sum(nsets(slist)),nlocs);
ii = 0;
sspoints = zeros(1,length(slist));

mood_all = cell(1,sum(nsets(slist)));
time_all = cell(1,sum(nsets(slist)));

freq = [4 8];
for sn = slist %[1,3,4,6,7,8,9,11,14,15,16] %already ran n.6   %Subjects showing enough variation in mood
    
    sub = subn(sn,:);
    [time,mood] = plot_mood(sub,false,true);
    close(gcf)
%     for k = 1:length(time)
%         F = griddedInterpolant(time{k},mood{k},'pchip');
%         xi = time{k}(1):time{k}(end);
%         mood{k} = F(xi);
%         time{k} = xi;
%     end
   
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    cd(data_path)
    name_list = dir;
    data_names = cell(1);
    
    time_all(ii+(1:length(time))) = time; 
    mood_all(ii+(1:length(time))) = mood;
    
    for k = 1:length(name_list)
        if strncmp(name_list(k).name, data_name, 18) && ~strcmp(name_list(k).name, '24201MMI_mmi3_proc1.ds')
            ii = ii+1;
            load([data_path,name_list(k).name,'/beamforming/VE_ICA_envelopes_',...
                num2str(freq(1)),'-',num2str(freq(2)),'Hz.mat'])
%             delete([data_path,name_list(k).name,'/beamforming/VE_ICA_enveopes.mat'])
            VE(ii,:) = VEhilb;         
            clear VEhilb
            
        end
    end

end


%% Concatenate data

VEall = cell(1,nlocs);
VEz = cell(1,nlocs);
posarray = false(1,nlocs); 
for ii = 1:length(sourcemodel.inside)
    a = false(1,size(VE,1));
    for jj = 1:size(VE,1)
        a(jj) = size(VE{jj,ii},2) > 1;
    end
    if all(a)
        VEall{ii} = cell2mat(VE(:,ii)');
        VEtemp = cell(1,size(VE,1));
        for k = 1:size(VE,1)            
            VEtemp{k} = zscore(VE{k,ii}); 
        end
        VEz{ii} = cell2mat(VEtemp)*std(VEall{ii})+mean(VEall{ii}); % normalize over subjects while maintaing variance and mean of signal in that voxel
        posarray(ii)  = true;    
    end
end

VEall = cell2mat(VEall')'; % location x time points
VEz =cell2mat(VEz')';
%% Z-score per subject
spoints = zeros(1,length(slist));
ii = find(posarray);
for jj = 1:size(VE,1)
    spoints(jj) = size(VE{jj,ii(1)},2);
end
% 
% spointss = [];
% a = cumsum(spoints);
% for ii = 1:length(spoints)
%     if spoints(ii)<500
%         if ii ==1
%             spointss = cat(2,spointss,1:a(ii));
%         else
%             spointss = cat(2,spointss,(a(ii-1)+1):a(ii));
%         end
%     end
% end
% VE = VEall;
% VE(spointss,:) = [];

% Z-score per subject
% a = nsets(slist); % number of sets per subject (included subs only)
% jj =0;
% spointss = zeros(1,length(slist));
% for ii = 1:length(slist)
%     spointss(ii) = sum( spoints(jj+ (1:a(ii)) ) );
%     jj = jj+a(ii);
% end
% 
% VEallz = VEall;
% a = cumsum(spointss);
% for ii = 1:length(slist)
%     if ii ==1
%         spoints = 1:a(ii);
%     else
%         spoints = (a(ii-1)+1):a(ii) ; 
%     end
%     VEallz(spoints,:) = zscore(VEallz(spoints,:)) ; 
% end
% VE = VEallz;
%%
[coeff, score] = pca(VEz);

% keep first 25 PCs
VEwhite = [];
VEwhite.time{1} = 1:size(VEall,1); 
VEwhite.trial{1} =  (score(:,1:30)*coeff(:,1:30)')';
% VEwhite.trial{1} =  VEz';
VEwhite.fsample = 1;
VEwhite.sampleinfo = [1 size(VEall,1)];
VEwhite.label = cell(size(VEall,2),1);
for ii = 1:size(VEall,2) % create unique channel names
VEwhite.label{ii} = ['VE',num2str(ii)];
end

cfg =[];
cfg.method = 'icasso';
cfg.icasso.method = 'fastica';
cfg.icasso.lastEig = 30; 
cfg.numcomponent = 100;
% cfg.fastica.numOfIC = 25;
cfg.channel = 'all';
cfg.icasso.Niter = 20;
comp = ft_componentanalysis(cfg, VEwhite);

% cfg =[];
% cfg.method = 'fastica';
% cfg.fastica.numOfIC = 25;
% cfg.fastica.lastEig = 30; 
% cfg.channel = 'all';
% comp = ft_componentanalysis(cfg, VEwhite);


save(['/data/MBDU/MEG_MMI3/results/mmi_ICA_envelopes_',num2str(freq(1)),...
    '-',num2str(freq(2)),'Hz.mat'],'comp')

%% Loading ICA analysis
freq = [4 8];
% freq = [8 13];
% freq = [13 30];
% freq = [35 90];

load(['/data/MBDU/MEG_MMI3/results/mmi_ICA_envelopes_',num2str(freq(1)),...
    '-',num2str(freq(2)),'Hz.mat'])

%% GLM

time = cell2mat(time_all);
mood = cell2mat(mood_all);

a = cumsum(spoints);
a(2:(end+1)) = a;
a(1) = 0;
ii = 1;
mdl = [];
for sn = slist   
%     y = cell2mat(mood_all(ii:(ii+nsets(sn)-1)));    
%     X = comp.trial{1}(:,(a(ii)+1):a(ii+nsets(sn)));   
%     ii = ii+nsets(sn);
%     mdl{sn} = fitglm(X',y','linear','Distribution','gamma');      
    y = cell2mat(mood_all(ii:(ii+nsets(sn)-1)));    
    y = zscore(y);
    X = comp.trial{1}(:,(a(ii)+1):a(ii+nsets(sn)));   
%     X = VEz((a(ii)+1):a(ii+nsets(sn)),:)'; 
    ii = ii+nsets(sn);
    mdl{sn} = fitglm(X',y','linear','Distribution','normal'); 
end

beta = zeros(26,length(slist));
ii = 0;
for sn = slist
    ii = ii+1;
    beta(:,ii) = mdl{sn}.Coefficients.Estimate;      
end

h = zeros(1,25);
p = zeros(1,25);
for ii = 1:25
[h(ii),p(ii)] = ttest(beta(ii+1,:)); % check mean is not zero, first beta is intercept
end

find(h)
[pv,ind] = sort(p);
%% Plot components
for ii = find(h)
% plot Indipendent components
compplot = zeros(nlocs,1); 
compplot(posarray) = comp.topo(:,ii)*sign(mean(beta(ii+1,:)));

sourceant =[];
sourceant.pow = compplot;
sourceant.dim = sourcemodel.dim;
sourceant.inside = sourcemodel.inside;
sourceant.pos = sourcemodel.pos;
cfg = [];
cfg.parameter = 'pow';
sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , mri);


crang = [-1 1]*0.02;
cfg = [];
cfg.method        = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourceant_Int);

set(gcf,'name',sprintf('Component %d, %d-%dHz, p=%.4f',ii,freq(1),freq(2),p(ii)))
% 
% figure(2); clf
% plot(comp.time{1}, comp.trial{1}(ii,:))
% xlim([a(14-1)+1 a(14)]); grid on
end

%%
% 
% cfg =[];
% cfg.method = 'pca';
% [coeff, score] = pca(VEhilb);
% score = comp_pca.trial{1}';
% compvar95 = cumsum(var(score,0,1))/sum(var(score,0,1)) <= 0.95;
% icomp = nnz(compvar95) ;
% clc
% fprintf('%d components for 95perc. of data variance\n',icomp)
% 
% if icomp>30
%     disp('Reducing ICA components to 30')
%     icomp = 30;
% end
% cfg =[];
% cfg.method = 'fastica';
% cfg.fastica.numOfIC = icomp;
% comp = ft_componentanalysis(cfg, data);
% 
% 
% figure
% cfg           = [];
% cfg.component = [1:icomp];       % specify the component(s) that should be plotted
% cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, comp)
% 
% 
% cfg          = [];
% cfg.channel  = [1:5]; % components to be plotted
% cfg.viewmode = 'component';
% cfg.layout   = 'CTF275.lay'; % specify the layout file that should be used for plotting
%     ft_databrowser(cfg, comp)


