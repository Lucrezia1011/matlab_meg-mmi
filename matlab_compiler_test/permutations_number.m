clc
cd /data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_old/


fit_parameters = {'mood';'E_LTA';'E_sum';'RPE_LTA';'RPE_sum'};

freqb = {'evoked';'delta';'theta';'alpha';'beta'};
eventb = {'cue';'choice';'outcome'};

for ii = 1:length(freqb)
    for jj = 1:length(eventb)
        for kk = 1:length(fit_parameters)
            
            if exist(sprintf('%s_%s/lmixmodel_%s/ROI_permute2.txt',freqb{ii},eventb{jj},fit_parameters{kk}),'file' )
                fprintf('%s_%s : %s\n',freqb{ii},eventb{jj},fit_parameters{kk});
                eval(sprintf('!wc -l %s_%s/lmixmodel_%s/ROI_permute2.txt ',freqb{ii},eventb{jj},fit_parameters{kk}));
            end
            
            if exist(sprintf('%s_%s/lmixmodel_%s/ROI_permute2/',freqb{ii},eventb{jj},fit_parameters{kk}),'dir' )
                fprintf('%s_%s : %s\n',freqb{ii},eventb{jj},fit_parameters{kk});
                eval(sprintf('!ls -1q %s_%s/lmixmodel_%s/ROI_permute2/ | wc -l',freqb{ii},eventb{jj},fit_parameters{kk}));
            end
        end
    end
end


%%


clc
cd /data/MBDU/MEG_MMI3/results/mmiTrial_aal_prep_mu5max/latent_vars_new/


fit_parameters = {'mood';'E_LTA';'E_sum';'RPE_LTA';'RPE_sum'};

freqb = {'evoked';'delta';'theta';'alpha';'beta'};
% eventb = {'cue';'choice';'outcome'};
eventb = {'outcome'};

for ii = 1:length(freqb)
    for jj = 1:length(eventb)
        for kk = 1:length(fit_parameters)
            
            if exist(sprintf('%s_%s/lme_%s/ROI_permute2.txt',freqb{ii},eventb{jj},fit_parameters{kk}),'file' )
                fprintf('%s_%s : %s\n',freqb{ii},eventb{jj},fit_parameters{kk});
                eval(sprintf('!wc -l %s_%s/lme_%s/ROI_permute2.txt ',freqb{ii},eventb{jj},fit_parameters{kk}));
            end
            
            if exist(sprintf('%s_%s/lme_%s/ROI_permute2/',freqb{ii},eventb{jj},fit_parameters{kk}),'dir' )
                fprintf('%s_%s : %s\n',freqb{ii},eventb{jj},fit_parameters{kk});
                eval(sprintf('!ls -1q %s_%s/lme_%s/ROI_permute2/ | wc -l',freqb{ii},eventb{jj},fit_parameters{kk}));
            end
        end
    end
end