function mmi_LTA_trials_sens_permute(latent_vars_name,n,npoints,fit_parameter,SD)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,npoints,fit_parameter,SD)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit
npoints = str2double(npoints);
N = str2double(n);

param_list = cell(N,1);
for nn = 1:N
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

load('neighbours.mat')
% datapath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/';

% meg_data_name = cat(2,data_path,meg_data_name);
% latent_vars_name = cat(2,data_path,latent_vars_name);

% cd /lscratch/$SLURM_JOBID

opts = detectImportOptions(latent_vars_name);
X = readtable(latent_vars_name,opts);

ntrials = max(X.trial)+1;% trials start from 0
Xp = X(:,1:2);
eval(sprintf('Xp.latentvar = X.%s;',fit_parameter));
clear X
ntot = size(Xp,1);

% Apos = cell(1,N);
% Aneg = cell(1,N);

% TFCE = zeros(N,npoints);
sub = unique(Xp.subject);

% Each subject is a block
% Seed for random number generator set in the swarm
SD = str2double(SD);
rng(SD)

rp = randperm(ntrials)-1; % trials start from 0
f = sign(rand(ntot,1)-0.5); % sign flip

%%
for nn = 1:N
    meg_data_name = ['meg_trials/sens_',param_list{nn},'.txt'];
%     megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 ntot-1])';
    megdata = dlmread(meg_data_name)';
    %     LME = table;
    LME = zeros(N,npoints);
    for z = 1:npoints
        X = Xp;
        X.MEG = megdata(:,z).*f;
        
        % Permute and sign flip MEG data within subject over all time
        for s = sub'
            trials = X.trial(X.subject == s); %find trial numbers per subject
            varx = X.MEG(X.subject == s); %find variables per subject
            trialn = ismember(rp',trials); % which trials are present
            [~,ind] = sort(rp(trialn)); %indeces of permuted trials numbers
            varxp = zeros(max(ind),1);
            varxp(ind) = varx;
            X.MEG(X.subject == s) = varxp;
        end
        
        lme = fitlme(X,'MEG ~ latentvar + (latentvar | subject) + (1|trial)');
        
        %     lme_names = lme.Coefficients.Properties.VarNames;
        
        %         lme_res = table;
        %         lme_res.ROI = nn;
        %         lme_res.index = z;
        %         lme_res.tStat = lme.Coefficients.tStat(2);
        %     for ii = 2:length(lme_names)
        %         eval(['lme_res.',lme_names{ii},' = lme.Coefficients.',lme_names{ii},'(2);'])
        %     end
        %         LME = cat(1,LME,lme_res);
        
        LME(nn,z) = lme.Coefficients.tStat(2);
    end
           
    fprintf('Done ROI %.0f\n',nn)
%     
end

E = 0.5; %0.5  % try and change the parameters
H = 2; %2
dh = 0.1;

tfce = matlab_tfce_transform_MEGtime(LME,H,E,dh,neighbours);
tfcen = matlab_tfce_transform_MEGtime(-LME,H,E,dh,neighbours);
tfce = tfce - tfcen;

%%
% Calculate peak maxima? If a cluster has bigger spatial extent than this
% then fine
Apos = max(tfce(:));
Aneg = min(tfce(:));

%
% outpath = [data_path,outpathn,'/lmixmodel_',fit_parameter,'/'];
% if ~exist(outpath,'dir')
%     mkdir(outpath)
% end

% outname = [outpath,'ROI_permute2.txt'];
% dlmwrite(outname,[Apos,Aneg],'-append')

dlmwrite('ROI_permute.txt',[Apos,Aneg])
% dlmwrite([meg_data_name,'_permute.txt'],LME)

end


function [tfced] = matlab_tfce_transform_MEGtime(T,H,E,dh,neighbours)
%MATLAB_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = matlab_tfce_transform(img,H,E,C,ndh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent,   H = 2
%   -- E extent exponent,   E = 0.5
%   -- A adjecency matrix of neighbouring channels
%   -- dh size of steps for cluster formation   dh = 0.1
% https://github.com/markallenthornton/MatlabTFCE/blob/master/matlab_tfce_transform.m

% set cluster thresholds
threshs = 0:dh:max(T(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(T(:));

% find connected components
vals = zeros(nvox,1);
if size(T,1)~=length(neighbours)
    T = T';
end
% cfg_neighb        = [];
% cfg_neighb.method = 'template';%'distance';
% cfg_neighb.template = 'CTF275_neighb.mat';
% cfg_neighb.channel = channels;
% neighbours        = ft_prepare_neighbours(cfg_neighb, hdr);

% channels = {neighbours.label};
% for n = 1:length(neighbours)
%     [~,~,iB] = intersect(neighbours(n).neighblabel, channels );
%     neighbours(n).neighbnum =iB;
% end

N = size(T,1);
for h = 1:ndh
    
      
    edgeA = cell(size(T,2),1);
    
    parfor tt = 1:size(T,2)-1
%      
        % Add spatial edges for next time points
        
        % Only need to check t+1 (t-1 is equivalently calculated in previous loop iteration)
%         Ttresh_m = T(:,tt-1)>=threshs(h);
        Ttresh = T(:,tt)>=threshs(h);
        Ttresh_p = T(:,tt+1)>=threshs(h);
        
        idx = find(Ttresh);
      
        neighboursT = neighbours(Ttresh);    
        edge_Add = cell(length(neighboursT),1);
        A = zeros(N);
        for n = 1:length(neighboursT)
            B = zeros(N,1);
            B(neighboursT(n).neighbnum) = 1;
            B(idx(n)) = 1;
            edge_p = find(Ttresh_p & B);
%             edge_m = find(Ttresh_m & B);
%             edge_Add{n,1} = [ones(nnz(edge_p)+nnz(edge_m),1)*idx(n) + N*(tt-1), ...
%                 [edge_p + N*(tt); edge_m + N*(tt-2)] ];
            edge_Add{n,1} = [ones(nnz(edge_p),1)*idx(n) + N*(tt-1), edge_p + N*(tt)];
            A(idx(n),Ttresh & B)  =1;                
        end
        Gt = graph(A);
        
        edgeA{tt} = [Gt.Edges.EndNodes + N*(tt-1); cell2mat(edge_Add)];      
        
    end
    
    tt = size(T,2);
    Ttresh = T(:,tt)>=threshs(h);
    
    idx = find(Ttresh);
    neighboursT = neighbours(Ttresh);
    A = zeros(N);
    for n = 1:length(neighboursT)
        B = zeros(N,1);
        B(neighboursT(n).neighbnum) = 1;
        B(idx(n)) = 1;
        A(idx(n),Ttresh & B)  =1;
    end
    Gt = graph(A);
    
    edgeA{tt} =Gt.Edges.EndNodes + N*(tt-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edgeA = cell2mat(edgeA);
    
        G = graph(edgeA(:,1),edgeA(:,2),1);
        nodes = unique(G.Edges.EndNodes);
        bins = conncomp(G);
        binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
        binnodes1 = [];
        n = 0; % only include with edges
        for binidx = 1:numel(binnodes)    
            if ismember(binnodes{binidx}, nodes) 
                n = n+1;
                binnodes1{n,1}  = binnodes{binidx};
            end
        end
    
        clustsize = zeros(nvox,1);

        voxpercc = cellfun(@numel,binnodes1);
        for c = 1:length(binnodes1)
            clustsize(binnodes1{c}) = voxpercc(c);
        end
        % calculate transform
        curvals = (clustsize.^E).*(threshs(h)^H);
        vals = vals + curvals;
    
    clc; fprintf('Calculating %.0f/%.0f\n',h,ndh)
end
tfced = NaN(size(T));
tfced(:) = vals.*dh;

end
