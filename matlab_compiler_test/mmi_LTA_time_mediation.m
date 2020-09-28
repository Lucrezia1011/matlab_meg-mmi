
% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit

%% Read Swarm output

nrois = 116;
data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_aal/';

% nrois = 269;
% data_path = '/data/MBDU/MEG_MMI3/results/mmiTrial_sensors_prep/';

param_list = cell(nrois,1);
for nn = 1:nrois
    n = num2str(nn);
    if size(n,2) == 1
        n = ['00',n];
    elseif size(n,2) == 2
        n = ['0',n];
    end
    param_list{nn} = n;
end

freq = 'evoked_outcome'; 
% freq = 'theta_outcome'; 

if strcmp(freq,'evoked_choice')
    npoints = 360;
    timew = [-.5,.7];
else
    npoints = 360;
    timew = [-.2,1];
end
meg_data_name = ['meg_trials_',freq,'.txt'];
latent_vars_name = ['latent_vars_',freq,'.csv'];


opts = detectImportOptions([data_path,latent_vars_name]);
Xv = readtable([data_path,latent_vars_name],opts);

fit_parameters = Xv.Properties.VariableNames(3:10);

%% Read data from file
Xfit = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = [];

    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s/lme_%s',data_path,freq,fit_parameter))
    
    for ii = 1:length(param_list)
        filename = sprintf('ROI_%s.csv',param_list{ii});
        opts = detectImportOptions(filename);
        x = readtable(filename,opts);
        x.index = x.index + (ii-1)*npoints; % adjust index
        X = cat(1,X,x);
    end
    Xfit{m} = X;
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement

clusteraal = cell(1,length(fit_parameters));
for m = 1:length(fit_parameters)
    X = Xfit{m};
%     if strncmp(freq,'evoked_',7)
        r = X.ROI;
%     else
%         [r,t] = ind2sub(size(meg),X.index); % find real ROI positions
%     end
    n = 0;
    clusteraal{m} = zeros(length(param_list),npoints);

    for iir = 1:length(param_list)

        x = X(r==iir,:);
        LME = x.tStat;
        
        TFCE = tfce2d(LME);

        clusteraal{m}(iir,:) = TFCE';
        % Equivalent to 
%         [tfced] = matlab_tfce_transform(LME,2,0.5,4,0.1);   
%         [tfcedn] = matlab_tfce_transform(-LME,2,0.5,4,0.1);   
%         tfced = tfced - tfcedn;
    end
    fprintf('Read parameter %d/%d\n',m,length(fit_parameters))
end


%% Calculate Threshold free cluster enhancement
close all
pv = 0.05;
plotopt = 't';

% Tested 3 predictors
M = 3;%length(fit_parameters);
fit_parameters_test = {'mood';'RPE_LTA';'RPE_sum'};
Xm = zeros(size(Xv,1),M);
for m = 1:M
    Xm(:,m) = Xv.(fit_parameters_test{m});
end
C = corr(Xm);
lambda = eig(C);
% Effective number of independent variables
Meff = 1 + (M-1)*(1 - var(lambda)/M);
alpha = 1 - (1 - pv)^(1/Meff);

% cannot combine over condition, because number of trials would be different  
% but can combine over predictor? 
ind = strfind(freq,'_');
freqb = freq(1:(ind-1));
clusternull = cell(length(fit_parameters),1);
condition = freq((ind+1):end);
for m = 1:length(fit_parameters)
    fit_parameter = fit_parameters{m};
    cd(sprintf('%s%s_%s/lme_%s',data_path,freqb,condition,fit_parameter))
    if exist('ROI_permute.txt','file')
        clusternull{m} = dlmread('ROI_permute.txt');      
    end
    if exist('ROI_permute','dir')
        nullnames = dir('ROI_permute');
        nullnames(1:2)=[];
        clusternull2 = cell(length(nullnames),1);
        for n = 1:length(nullnames)
            if nullnames(n).bytes == 0
                delete(['ROI_permute/',nullnames(n).name])
            else
                clusternull2{n} = dlmread(['ROI_permute/',nullnames(n).name]);
            end
        end
        clusternull{m} = cat(1,clusternull{m},cell2mat(clusternull2));
    end  
end

clusternull = cell2mat(clusternull);
if isempty(clusternull)
    clusternull = cell2mat(clusteraal);
    alpha = 0.05/116;  
end

snull = sort(abs(clusternull(:)));
Nnull = size(snull,1);

uplim =  snull(ceil(Nnull*(1-alpha)),1);
lowlim = -uplim;

%%
aal_labels = readcell('~/labels_AAL116_MNIv4.csv');
time = linspace(-.2,1,npoints);
Xv.RPE_abs = abs(Xv.RPE);


% mmi_LTA_aal_trials(meg_data_name,latent_vars_name,n,fit_parameter)
% meg_data_name: name of text file with meg data (nchans*npoints, ntrials)
% latent_vars_name:  name of table file with latenet varibles per trial
% n = meg row to fit (starts from 0) total npoints*nrois = 300*116 =0:34799
% fit_parameter = name of latent variable to fit

fit_parameter1 = 'RPE_sum';
fit_parameter2 = 'E';%'RPE_abs';%'E';

ns = zeros(1,nrois);
for n  =1:nrois
    ns(n) = any(abs(clusteraal{strcmp(fit_parameters,fit_parameter1)}(n,:))>uplim);
end
ns = find(ns);


subs = unique(Xv.subject)';
for s = subs
    Xv.(fit_parameter1)(Xv.subject==s) = zscore(Xv.(fit_parameter1)(Xv.subject==s));
    Xv.(fit_parameter2)(Xv.subject==s) = zscore(Xv.(fit_parameter2)(Xv.subject==s));
end

%%
% Mediation
lme2 = fitlme(Xv,sprintf('%s ~ %s +(%s | subject) + (1|trial)',...
    fit_parameter2,fit_parameter1,fit_parameter1));

mediation = cell(1,length(ns));
for ni = 1:length(ns)
nn = ns(ni);


zz = find(abs(clusteraal{strcmp(fit_parameters,fit_parameter1)}(nn,:))>uplim);


megdata = dlmread([data_path,meg_data_name],',',[(nn-1)*npoints 0 nn*npoints-1 size(Xv,1)-1])';

% Test that there is no complete mediation
M = cell(1,length(zz));
parfor z = 1:length(zz)
    ltv = Xv;
    ltv.MEG = megdata(:,zz(z));
    for s = subs
        % z-score over subject
        ltv.MEG(ltv.subject==s) = zscore(ltv.MEG(ltv.subject==s));
    end
    
    % Main effect
    lme1 = fitlme(ltv,sprintf('MEG ~ %s +(%s | subject) + (1|trial)',...
        fit_parameter1,fit_parameter1));
%     % Mediation
%     lme2 = fitlme(ltv,sprintf('%s ~ %s +(%s | subject) + (1|trial)',...
%         fit_parameter2,fit_parameter1,fit_parameter1));
    
    % combined
    lme3 = fitlme(ltv,sprintf('MEG ~ %s + %s + (%s + %s | subject) + (1|trial)',...
        fit_parameter1,fit_parameter2,fit_parameter1,fit_parameter2));
    
    LME = [];
    % Step 1
    LME.xy_tStat = lme1.Coefficients.tStat(2); % Equivalent of Estimate/SE
    LME.xy_pValue = lme1.Coefficients.pValue(2);
    
    % Step 2
    LME.xm_tStat = lme2.Coefficients.tStat(2);
    LME.xm_pValue =lme2.Coefficients.pValue(2);
    
    % Step 3
%     LME.xym_Name = lme3.Coefficients.Name;
%     LME.xym_tStat = lme3.Coefficients.tStat(3);
%     LME.xym_pValue = lme3.Coefficients.pValue(3);
    
    LME.c = lme1.Coefficients.Estimate(2);
    LME.a = lme2.Coefficients.Estimate(2);
    LME.c1 = lme3.Coefficients.Estimate( strcmp(lme3.Coefficients.Name,fit_parameter1));
    LME.b = lme3.Coefficients.Estimate( strcmp(lme3.Coefficients.Name,fit_parameter2));
   
    
    % The value of the mediated or indirect effect estimated by taking the 
    % difference in the coefficients, ĉ − ĉ′, corresponds to the reduction 
    % in the independent variable effect on the dependent variable when 
    % adjusted for the mediator. To test for significance, the difference 
    % is then divided by the standard error of the difference and the ratio 
    % is compared to a standard normal distribution.
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819368/
    
    % Step 4
    LME.mediation = (LME.c-LME.c1)/sqrt(lme1.Coefficients.SE(2)^2 + lme3.Coefficients.SE( strcmp(lme3.Coefficients.Name,fit_parameter1))^2);
    LME.ab = LME.a*LME.b;
   
    % Sobel Test:  SE_{ab} = b^2*SEb^2 + a^2*SEa^2
    LME.SE_ab = sqrt( LME.a^2 * lme2.Coefficients.SE(2)^2 + ...
        LME.b^2 * lme3.Coefficients.SE( strcmp(lme3.Coefficients.Name,fit_parameter2))^2 );
    
   
    %   The amount of mediation is called the indirect effect.   Note that the
    %   total effect = direct effect + indirect effect; 
    %   or using symbols:  c = c' + ab
    %  the two are only approximately equal for multilevel models, logistic 
    %  analysis and structural equation modeling with latent variables.  
    %  For such models, it is probably inadvisable to compute c from Step 1, 
    %  but rather c or the total effect should be inferred to be c' + ab and 
    %  not directly computed.


    % Joint Significance of Paths a and b: Test that a and b are not zero.
    % i.e. xm_tStat & lme3.Coefficients.tStat( strcmp(lme3.Coefficients.Name,'E'))
    LME.jointSig = LME.xm_pValue<0.05 & ...
        lme3.Coefficients.pValue( strcmp(lme3.Coefficients.Name,fit_parameter2))<0.05;
    % Z-test > 1.96 in absolute value is significant at the .05 level
    LME.abz = LME.ab/LME.SE_ab;    
    M{z} = LME;
    clc
    fprintf('Done %.1f\n',z/length(zz)*100)
    
    
end
mediation{ni} = cell2mat(M);

end
%Most contemporary analysts believe that the essential steps in establishing 
% mediation are Steps 2 and 3. Certainly, Step 4 does not have to be met 
% unless the expectation is for complete mediation.  In the opinion of most 
% though not all analysts, Step 1 is not required.

% Monte Carlo http://quantpsy.org/medmc/medmc.htm

%%

for ni = 1:length(ns)
nn = ns(ni);
figure(ni); clf; 
set(gcf,'name',aal_labels{nn},'color','w','position',[566   268   415   731])

T2 = Xfit{strcmp(fit_parameters,fit_parameter2)}.tStat(Xfit{1}.ROI== nn);
T1 = Xfit{strcmp(fit_parameters,fit_parameter1)}.tStat(Xfit{1}.ROI== nn);
subplot(311); plot(time,T1)
hold on
plot(time,T2)
title([aal_labels{nn},' Fixed effect'])
grid on
ylabel('t-value'); xlabel('time (s)')
legend(fit_parameter1,fit_parameter2,'location','best')

subplot(312)
plot(time,clusteraal{strcmp(fit_parameters,fit_parameter1)}(nn,:))
hold on
plot(time,clusteraal{strcmp(fit_parameters,fit_parameter2)}(nn,:))
plot([-.2,1],[1 1]*uplim,'k')
plot([-.2,1],[1 1]*lowlim,'k')
title(['Fixed effect (temporal clustering)'])
legend(fit_parameter1,fit_parameter2,'p<0.05','location','best')
grid on
zz = find(abs(clusteraal{strcmp(fit_parameters,fit_parameter1)}(nn,:))>uplim);
ylabel('TFCE'); xlabel('time (s)')

subplot(313)
abz = zeros(1,length(zz));
m = zeros(1,length(zz));
abc = zeros(1,length(zz));
for z  =1:length(zz)
    abz(z) = mediation{ni}(z).abz;
    m(z) = mediation{ni}(z).mediation;
    abc(z) = mediation{ni}(z).ab/mediation{ni}(z).c;
end
plot(1:length(zz),abz,'-x','MarkerSize',12)
hold on
plot(1:length(zz),m,'-x','MarkerSize',12)
plot(1:length(zz),abc,'-x','MarkerSize',12)
plot([0,length(zz)+1],[1 1]*1.96,'k')
plot([0,length(zz)+1],[1 1]*-1.96,'k')
if nnz(zz)>8
    set(gca,'XTick',1:5:nnz(zz),'XTickLabel',round(time(zz(1:5:nnz(zz)))*1000))
else
    set(gca,'XTick',1:nnz(zz),'XTickLabel',round(time(zz)*1000))
end
legend('a*b','c-c''','ab/c','location','best')
ylabel('z-value'); xlabel('time (ms)')
title(sprintf('%s --> %s --> MEG',fit_parameter1,fit_parameter2))
grid on
saveas(gcf,sprintf('~/matlab/figures/%s_mediator-%s_%s.png',...
    fit_parameter1,fit_parameter2,aal_labels{nn}))
end


