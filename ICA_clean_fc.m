clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

for ii = [15,16] %[1,3,4,6,7,8,9,11,14,15,16] % Subjects showing enough variation in mood
    sub = subn(ii,:);
    
    data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
    cd(data_path)
    data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 0.5-300 Hz to adjust baseline
    
    if ~exist(data_name,'dir')
        data_name = [sub,'MMI_mmi3_proc1.ds'];
    end
    
    if strcmp(sub, '24201') || strcmp(sub, '22695')
        data_name = [sub,'MMI_mmi3_proc2.ds'];
    end
  
    %% Clean Data
    
    % Unsure of best pre-processing for baseline in sub 24071
    cfg = [];
    cfg.dataset = data_name;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';
    % cfg.demean = 'yes';
    % cfg.bpfilter = 'yes';
    % cfg.bpfreq = [1 150];
    data = ft_preprocessing(cfg);
    f = data.fsample;
    
    if exist([data_path,'ICA_artifacts.mat'],'file')
        load([data_path,'ICA_artifacts.mat']);
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
        
        if icomp>40
            disp('Reducing ICA components to 40')
            icomp = 40;
        end
        cfg =[];
        cfg.method = 'fastica';
        cfg.fastica.numOfIC = icomp;
        comp = ft_componentanalysis(cfg, data);
        
        
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
        
        try
            figure;
            plot(abs(corr(eog',comp.trial{1}')))
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
%     
%     cfg           = [];
%     cfg.component = 1:length(comps.label);
%     data_clean    = ft_rejectcomponent(cfg, comps,data);
    
    
   
end