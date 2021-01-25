function mmi_LTA_powergrid_permute(meg_data_name,latent_vars_name,fit_parameter,SD)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit


% datapath = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/';

% meg_data_name = cat(2,data_path,meg_data_name);
% latent_vars_name = cat(2,data_path,latent_vars_name);

opts = detectImportOptions(latent_vars_name);
X = readtable(latent_vars_name,opts);
% gridin = dlmread('mni_grid.txt');

fid=fopen(meg_data_name);
g = textscan(fid,'%s','delimiter','\n');
fclose(fid);
ntot = length(g{1});
clear g
% if (nn*npoints)<ntot
%     megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 size(X,1)-1])';
% else
%     megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 ntot-1 size(X,1)-1])';
%     npoints = size(megdata,2);
% end

% read all data
megdata = dlmread(meg_data_name,',',[0 0 ntot-1 size(X,1)-1])';
npoints = size(megdata,2);

% Includes intercep random effect for subject
% Uncorrelated random effect for intercept and mood grouped by subject,
% plus random effect for intercept given trial

lme_formula = sprintf('MEG ~ %s +(%s|subject) + (1|trial) + (1|recording)',fit_parameter,fit_parameter);

% Each subject is a block
% Seed for random number generator set in the swarm
SD = str2double(SD);
rng(SD)
trials = unique(X.trial);
rp = randperm(length(trials));
rp = trials(rp);

LME = cell(1,npoints);
for z = 1:npoints
    X.MEG = megdata(:,z);
    
    % Permute MEG data within recording over all time
    for s = 1:max(X.recording)
        trials = X.trial(X.recording == s); %find trial numbers per subject
        varx = X.MEG(X.recording == s); %find variables per subject
        trialn = ismember(rp,trials); % which trials are present
        [~,ind] = sort(rp(trialn)); %indeces of permuted trials numbers
        varxp = zeros(max(ind),1);
        varxp(ind) = varx;
        X.MEG(X.recording == s) = varxp;
    end
    
    lme = fitlme(X,lme_formula); %Fixed effects for RPE

%     lme.Coefficients.Properties.VarNames;
%     
%     lme_names = lme.Coefficients.Properties.VarNames;
%      
%     lme_res = table;
%     lme_res.index = z;
%     for ii = 2:length(lme_names)
%         lme_res.(lme_names{ii}) = lme.Coefficients.(lme_names{ii})(2);
%     end
%     LME = cat(1,LME,lme_res);
    
    LME{z} = lme.Coefficients.tStat(2);
    clc
%     fprintf('LME done %.1f\n',z/npoints*100);
    
end
% 
% img= zeros(size(gridin));
% img(gridin==1) = cell2mat(LME);
% % sourcemodel.dim
% img = reshape(img,[32 39 34]);
% tfce = tfce_transform(img,2,0.5,26,0.1);%one tailed,positive only
% tfce = tfce(:);
% tfce = tfce(gridin==1);


% dlmwrite('grid_permute.txt',tfce); % max(tfce)== largest cluster; keep all
dlmwrite('grid_permute.txt',cell2mat(LME)); % save without tfce transform

end


function [tfced] = tfce_transform(img,H,E,C,dh)
%MATLAB_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = matlab_tfce_transform(img,H,E,C,ndh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent,   H = 2
%   -- E extent exponent,   E = 0.5
%   -- C connectivity,      C = 26
%   -- dh size of steps for cluster formation   dh = 0.1
% https://github.com/markallenthornton/MatlabTFCE/blob/master/matlab_tfce_transform.m

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(img(:));

% find connected components
vals = zeros(nvox,1);
cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,img,x),C), threshs);
for h = 1:ndh
    clustsize = zeros(nvox,1);
    ccc = cc(h);
    voxpercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:ccc.NumObjects
        clustsize(ccc.PixelIdxList{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(img));
tfced(:) = vals.*dh;

end