clear all
close all
clc


subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

%% Find common channels
sn = 2;
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds'];
h = ft_read_sens(data_name,'senstype','MEG');
label2 = h.label(strncmp(h.label,'M',1));

sn = 4;
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)
data_name = [sub,'MMI_mmi3_proc.ds'];
h = ft_read_sens(data_name,'senstype','MEG');
label4 = h.label(strncmp(h.label,'M',1));

channels = intersect(label2,label4);

%%

param_list = [];

zz= 0;
for sn =[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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
        zz = zz +1;
        param_list{zz} = data_names{runs};  
    end
end

% Run for all channels in common in order to run lmixm indipendently for
% all sensors
for ii = 1:length(param_list)
    mmi_LTA_sensors_prep(param_list{ii},channels)
end

%%

% save(save_name ,'Ybeta','Ytheta','Y','ltvall')

Yall = [];
ltv = [];
% Ym = [];
YTheta = [];
YBeta = [];
ntrials = [];
for sn =[1:7,11,14:16] %[1:6,8,9,14] % all subjects with continuos recordings and latent variables
        
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
        data_name = data_names{runs};
        n = str2double(data_name(end-3));
        if isnan(n)
            n = 1;
        end
        
        load(['/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/',sub,'_',num2str(n)])
        if isempty(ltv)
            ltv = ltvall;
        else
            ltv(end+(1:size(ltvall,1)),:) = ltvall;
        end
        % Evoked responses 200Hz sampling
        Yall = cat(3,Yall,Y);
        
        % Induced responses 100Hz sampling
        YTheta = cat(3,YTheta,Ytheta);
        YBeta = cat(3,YBeta,Ybeta);
        % Check sign of evoked responses
%         Ym = cat(3,Ym,mean(Y,3));
        ntrials = cat(1,ntrials,size(Y,3)); 
    end
end

ltv = flipud(ltv);

%% Sign flip for source localized evoked responses
% S = zeros(size(Yall));
% npoints = size(Ym,2);
% for ii = 1:size(Ym,1)
%     Yu = squeeze(Ym(ii,:,:))';
% %     [Yu,IA,IC] = unique(Y','rows');
%     Yu = zscore(Yu,0,2);
%    
%     %    % first pass
%     Yav = [];
%     C = corr(Yu');
%     C = triu(C,1);
%     C( abs(C-1) < 0.01) = 0;
%     [~,ind] = sort(abs(C(:)),'descend');
%     m = sign(C(ind(1)));
%     [i,j] = ind2sub(size(C),ind(1));
%     Yav = cat(1,Yav, Yu(i,:), m*Yu(j,:));
%     s = zeros(size(C,1),1);
%     s(i) = 1;
%     s(j) = m;
%     
%     while nnz(s) < size(Yu,1)
%         C = corr(mean(Yav,1)',Yu');
%         [~,ind] = sort(abs(C(:)),'descend');
%         inds = find(s);
%         z = 1;
%         while any(ind(z) == inds)
%             z = z+1;
%         end
%         
%         m = sign(C(ind(z)));
%         s(ind(z)) = m;
%         Yav =  cat(1,Yav, m*Yu(ind(z),:));
%     end
%     
% %     time = linspace(-0.5,1,size(Yav,2));
% %     figure; plot(time,mean(Yav,1))
%   
%     zz = 0;
%     for z = 1:size(Yu,1)
%         S(ii,:,zz+(1:ntrials(z))) = s(z);
%         zz = zz + ntrials(z);
%     end
%     
% end
% 
% Yall = Yall.*S;
Yall = reshape(Yall,size(Yall,1)*size(Yall,2),size(Yall,3)); % nrois * npoints * ntrials
YBeta = reshape(YBeta,size(YBeta,1)*size(YBeta,2),size(YBeta,3)); % nrois * npoints * ntrials
YTheta = reshape(YTheta,size(YTheta,1)*size(YTheta,2),size(YTheta,3)); % nrois * npoints * ntrials

%% Write data
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/meg_trials.txt',Yall);
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/meg_trials_beta.txt',YBeta);
dlmwrite('/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/meg_trials_theta.txt',YTheta);


writetable(ltv,'/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/latent_vars.csv');

