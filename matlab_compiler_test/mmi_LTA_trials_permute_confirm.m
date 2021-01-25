function mmi_LTA_trials_permute_confirm(meg_data_name,latent_vars_name,n,npoints,fit_parameter,SD)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,npoints,fit_parameter,SD)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit
npoints = str2double(npoints);
N = str2double(n); 

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

nperms = 20;

TFCE = zeros(nperms,npoints);
sub = unique(Xp.subject);

% Each subject is a block
% Seed for random number generator set in the swarm
SD = str2double(SD);
rng(SD)

nperms = 10;
%%

megdata = dlmread(meg_data_name,',',[(N-1)*npoints 0 N*npoints-1 ntot-1])';
for nn = 1:nperms
    
    % new randomization
    rp = randperm(ntrials)-1; % trials start from 0
%     f = sign(rand(ntot,1)-0.5); % sign flip
    
    X = Xp;
    megshuffle = megdata; % megdata.*f;  % no sign flipping
    % Permute and sign flip MEG data within subject over all time
    for s = sub'
        trials = X.trial(X.subject == s); %find trial numbers per subject
        varx = megshuffle(X.subject == s,:); %find variables per subject
        trialn = ismember(rp',trials); % which trials are present
        [~,ind] = sort(rp(trialn)); %indeces of permuted trials numbers
        varxp = zeros(max(ind),npoints);
        varxp(ind,:) = varx;
        megshuffle(X.subject == s,:) = varxp;
    end
    
    %     LME = table;
    LME = cell(1,npoints);
    for z = 1:npoints
        
        X.MEG = megshuffle(:,z);
        
        
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
        
        LME{z} = lme.Coefficients.tStat(2);
    end
    
    LME = cell2mat(LME)';
       
    
    E = 0.5;
    H = 2;
    S=regionprops(LME>0,'PixelIdxList','PixelList');
    TFCEp = tfextent(LME,S,E,H);
    
    S=regionprops(LME<0,'PixelIdxList','PixelList');
    TFCEn = tfextent(-LME,S,E,H);
    
    TFCEp(TFCEn>0) = -TFCEn(TFCEn>0);
    
    TFCE(nn,:) = TFCEp;
    clc
    fprintf('Done permutation %.0f\n',nn)
    
end
%%
% Calculate peak maxima? If a cluster has bigger spatial extent than this
% then fine
Apos = max(TFCE,[],2);
Aneg = min(TFCE,[],2);

%
% outpath = [data_path,outpathn,'/lmixmodel_',fit_parameter,'/'];
% if ~exist(outpath,'dir')
%     mkdir(outpath)
% end

% outname = [outpath,'ROI_permute2.txt'];
% dlmwrite(outname,[Apos,Aneg],'-append')

dlmwrite('ROI_permute.txt',[Apos,Aneg])

end


function TFCE = tfextent(Tst,S,E,H)

TFCE = zeros(size(Tst));
% convert into area
dx = 1;

for ii = 1:length(S)
    a = S(ii).PixelIdxList;
    
    h = Tst(a);
    
    dh = 0.1;
    for p = 1:length(h)
        A = 0;
        hp = dh;
        while hp < h(p)
            e = hp <= h;
            de = diff(e);
            de(end+1,:) = 0;
            e = e & ~de;
            for t = 1:length(e)
                if e(t) == 0 && t<p
                    e(1:t) = 0;
                elseif e(t) == 0 && t>p
                    e(t:end) = 0;
                end
            end
            % e =  time extent
            
            eE = (sum(e)*dx)^E;
            A = A + eE*(hp^H)*dh;
            hp = hp + dh;
        end
        TFCE(a(p)) = A;
    end
end
end