function mmi_LTA_trials_permute2(meg_data_name,latent_vars_name,n,npoints,fit_parameter,SD)
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
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

TFCE = zeros(N,npoints);
sub = unique(Xp.subject);

% Each subject is a block
% Seed for random number generator set in the swarm
SD = str2double(SD);
rng(SD)

rp = randperm(ntrials)-1; % trials start from 0
f = sign(rand(ntot,1)-0.5); % sign flip

%%
for nn = 1:N
    
    megdata = dlmread(meg_data_name,',',[(nn-1)*npoints 0 nn*npoints-1 ntot-1])';
  
%     LME = table;
    LME = cell(1,npoints);
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
        
        LME{z} = lme.Coefficients.tStat(2);        
    end
    
    LME = cell2mat(LME);
%     F = griddedInterpolant(1:npoints,LME,'pchip');
    dx = 1;
%     samples = 1:dx:npoints;
%     LMEi = F(samples);
 
    indp = LME>0;
    ii = find(diff(indp));
    if indp(1) == 0 && indp(end) == 0
        a = cell(1,length(ii)/2);
        for jj = 1:length(ii)/2
            a{jj} = ii(jj*2-1)+1:ii(jj*2);
        end
    elseif indp(1) == 1 && indp(end) == 0
        a = cell(1,ceil(length(ii)/2));
        a{1} = 1:ii(1);
        for jj = 1:floor(length(ii)/2)
            a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
        end
    elseif indp(1) == 0 && indp(end) == 1
        a = cell(1,ceil(length(ii)/2));
        for jj = 1:floor(length(ii)/2)
            a{jj} = ii(jj*2-1)+1:ii(jj*2);
        end
        a{end} = ii(end)+1:npoints;
    else
        a = cell(1,length(ii)/2+1);
        a{1} = 1:ii(1);
        for jj = 1:floor(length(ii)/2)-1
            a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
        end
        a{end} = ii(end)+1:npoints;
    end
    
    % convert into area
    for ii = 1:length(a)
        h = LME(a{ii});
        
        dh = 0.1;
        for p = 1:length(h)
            A = 0;
            hp = dh;
            while hp < h(p)
                e = hp <= h;
                de = diff(e);
                de(end+1) = 0;
                e = e & ~de;
                for t = 1:length(e)
                    if e(t) == 0 && t<p
                        e(1:t) = 0;
                    elseif e(t) == 0 && t>p
                        e(t:end) = 0;
                    end
                end
                eE = (sum(e)*dx)^0.5;
                
                A = A + eE*(hp^2)*dh;
                hp = hp + dh;
            end
            TFCE(nn,p+a{ii}(1)-1) = A;
        end   
    end
      
    indp = LME<0;
    ii = find(diff(indp));
    if indp(1) == 0 && indp(end) == 0
        a = cell(1,length(ii)/2);
        for jj = 1:length(ii)/2
            a{jj} = ii(jj*2-1)+1:ii(jj*2);
        end
    elseif indp(1) == 1 && indp(end) == 0
        a = cell(1,ceil(length(ii)/2));
        a{1} = 1:ii(1);
        for jj = 1:floor(length(ii)/2)
            a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
        end
    elseif indp(1) == 0 && indp(end) == 1
        a = cell(1,ceil(length(ii)/2));
        for jj = 1:floor(length(ii)/2)
            a{jj} = ii(jj*2-1)+1:ii(jj*2);
        end
        a{end} = ii(end)+1:npoints;
    else
        a = cell(1,length(ii)/2+1);
        a{1} = 1:ii(1);
        for jj = 1:floor(length(ii)/2)-1
            a{jj+1} = ii(jj*2)+1:ii(jj*2+1);
        end
        a{end} = ii(end)+1:npoints;
    end
    % convert into area
    for ii = 1:length(a)
        h = -LME(a{ii});
        
        dh = 0.1;
        for p = 1:length(h)
            A = 0;
            hp = dh;
            while hp <= h(p)
                e =  hp <= (h+dh/2);
                de = diff(e);
                de(end+1) = 0;
                e = e & ~de;
                for t = 1:length(e)
                    if e(t) == 0 && t<p
                        e(1:t) = 0;
                    elseif e(t) == 0 && t>p
                        e(t+1:end) = 0;
                    end
                end
                
                %                 hold on
                %                 plot(find(e),ones(1,nnz(e))*hp)
                eE = (sum(e)*dx)^0.5;
                
                A = A + eE*(hp^2)*dh;
                hp = hp + dh;
            end
            TFCE(nn,p+a{ii}(1)-1) = -A;
        end
    end
    clc
    fprintf('Done ROI %.0f\n',nn)
end
 %%
% surf(linspace(-.2,1,npoints),1:N,TFCE)
% Calculate peak maxima? If a cluster has bigger spatial extent than this
% then fine
Apos = max(TFCE(:));
Aneg = min(TFCE(:));

% 
% outpath = [data_path,outpathn,'/lmixmodel_',fit_parameter,'/'];
% if ~exist(outpath,'dir')
%     mkdir(outpath)
% end

% outname = [outpath,'ROI_permute2.txt'];
% dlmwrite(outname,[Apos,Aneg],'-append')

dlmwrite('ROI_permute2.txt',[Apos,Aneg])

end