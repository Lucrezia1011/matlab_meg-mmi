clear all
% close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults
addpath('~/fieldtrip-20190812/fieldtrip_private')


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

sn = 6; %[1,3,4,6,7,8,9,11,14,15,16] %Subjects showing enough variation in mood
sub = subn(sn,:);
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

dataprep = mmi_dataprep(sub);

%%
data = dataprep.data;
filt_order = []; % default
% data_filt = ft_preproc_bandpassfilter(data_clean.trial{1}, data.fsample, [0 35],filt_order,'but')';
data_filt = ft_preproc_lowpassfilter(data.trial{1}, data.fsample, 35 ); 

data.trial{1} = data_filt;
clear data_filt

%% Read events

bv_match = match_triggers_fc(dataprep.dataname);

answer_match = bv_match.answer;
choice_match =  bv_match.choice;
outcome_match  = bv_match.outcome;
mood_match = bv_match.ratemood;
blockmood_match = bv_match.blockmood;
slider_match = bv_match.slider;
blockslider_match = bv_match.blockslider;
ITI_match = bv_match.ITI ;
buttonpress = bv_match.buttonpress;

pRPE = outcome_match.win == 1 ;
nRPE = outcome_match.win == -1 ;
pRPE_sample = outcome_match.sample(pRPE);
nRPE_sample = outcome_match.sample(nRPE);

%% Minimum norm

% cfg.covariance = 'yes';
% cfg.covariancewindow = [-inf 0];
% avgpress = ft_timelockanalysis(cfg, datab);
% 
% cfg               = [];
% cfg.method        = 'mne';
% cfg.sourcemodel   = dataprep.leadfield;
% cfg.headmodel     = dataprep.headmodel;
% cfg.latency       = [0.035];
% cfg.mne.prewhiten = 'yes';
% cfg.mne.lambda    = 3;
% cfg.mne.scalesourcecov = 'yes';
% 
% sourceFC          = ft_sourceanalysis(cfg,avgpress);
% 
% 
% sourceant =[];
% sourceant.pow = mean(sourceFC.avg.pow,2);
% sourceant.dim = sourceFC.dim;
% sourceant.inside = sourceFC.inside;
% sourceant.pos = sourceFC.pos;
% cfg = [];
% cfg.parameter = 'pow';
% sourceant_Int  = ft_sourceinterpolate(cfg, sourceant , dataprep.mri);
% 
% 
% crang = [];
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter = 'pow';
% cfg.maskparameter = 'pow';
% cfg.funcolormap  = 'auto';
% cfg.funcolorlim   = crang;
% cfg.opacitylim = crang;
% ft_sourceplot(cfg, sourceant_Int);



%% Optimized beamforer, from Bill's paper
% Q(1:size(L,2)) = {0};
% P(1:size(L,2)) = {0};
% 
% % case 'singleshell'
% meansphereorigin = mean(dataprep.headmodel.bnd.pos,1);
% % the angles are the same for all dipole locations
% all_angles = linspace(0,pi,180); 
% 
% sdata = svd(data.trial{1}');
% noise = sdata(end-icacomps);
% parfor ii = 1:length(L)
%   optim_options = optimset('Display', 'final', 'TolX', 1e-3, 'Display', 'off');
%     if ~isempty(W{ii})
%                
%         lf = L{ii};
%         vox_pos = dataprep.leadfield.pos(ii,:);
%         
%                
%         % perform a non-linear search for the optimum orientation
%         [tanu, tanv] = calctangent(vox_pos - meansphereorigin); % get tangential components
%             
%         % get a decent starting guess
%         all_costfun_val = zeros(size(all_angles));
%         for i=1:length(all_angles)
%           costfun_val        = SAM_costfun(all_angles(i), vox_pos, tanu, tanv, lf, Cr, inv(Cr), noiseC);
%           all_costfun_val(i) = costfun_val;
%         end
%         [junk, min_ind] = min(all_costfun_val);
%  
%         
%         try
%             [opt_angle, fval, exitflag, output] = fminsearch(@SAM_erf_costfun, ...
%                 all_angles(min_ind), optim_options, vox_pos, tanu, tanv, lf, Cr, inv(Cr), noise, noiseC, dataerf);
%             MDip        = settang(opt_angle, tanu, tanv);
%             MagDip      = sqrt(dot(MDip,MDip));
%             opt_vox_or  = (MDip/MagDip)';
%         catch %if optimisaiton does not converge
%                   
% %             [v,d] = svd(lf'/Cr*lf);
% %             d = diag(d);
% %             if d(3) < 1
% %                 jj = 2; % The minumum singular value is degenerate
% %             else
% %                 jj = 3;
% %             end  
% %             opt_vox_or = v(:,jj);
%             
%             opt_angle = all_angles(min_ind);
%             MDip        = settang(opt_angle, tanu, tanv);
%             MagDip      = sqrt(dot(MDip,MDip));
%             opt_vox_or  = (MDip/MagDip)';
%         end
%         lfo = lf* opt_vox_or;      
%         w = Cr\lfo / (lfo'/Cr*lfo) ;
%         
%         wnorm = w/sqrt( sum( (w*noise(1,1)).^2) );
%         dataloc = zeros(ntrials,nsamples);
%         for tt = 1:ntrials
%             dataloc(tt,:) = wnorm'*dataerf.trial{tt};
%         end
%         P{ii} = mean(dataloc,1);
%         
%         
%     end 
%     
%     if mod(ii,300) == 0
%         fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
%     end
%     
% end

%% Beamfomer

C = cov(data.trial{1}');
E = svd(C); 
nchans = length(data.label);
icacomps = length(data.cfg.component);
noiseC = eye(nchans)*E(end-icacomps); % ICA eliminates from 2 to 4 components

Cr = C + noiseC; % need to normalise because of ICA

L = dataprep.leadfield.leadfield;
VEp(1:size(L,2)) = {0};
VEn(1:size(L,2)) = {0};
VEg(1:size(L,2)) = {0};
VEs(1:size(L,2)) = {0};

W(1:size(L,2)) = {0};

% dataerf = define_trials(buttonpress,data,bv_match.time,[-0.2 1]);

Rsamples = outcome_match.sample(outcome_match.win~=0);
dataerf = define_trials(Rsamples,data,bv_match.time,[-0.2 1]);


sdata = svd(cell2mat(dataerf.trial)');
noise = sdata(end-icacomps);


ntrials = length(dataerf.time);
nsamples = length(dataerf.time{1});

dataerfc = cell2mat(dataerf.trial);

parfor ii = 1:length(L)
    lf = L{ii}; % Unit 1Am
    if ~isempty(lf)
        
        % %  G O'Neill method, equivalent to ft
        [v,d] = svd(lf'/Cr*lf);
        d = diag(d);
        jj = 2;
%         if d(3) < 1
%             jj = 2; % The minumum singular value is degenerate
%         else
%             jj =3;
%         end
        lfo = lf*v(:,jj); % Lead field with selected orientation
        

        w = Cr\lfo / (lfo'/Cr*lfo) ;      
       
%         wnorm = Cr\lfo / sqrt(lfo'/(Cr^2)*lfo) ;
        
        W{ii} = w;
                
       
    end
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
    
end
clc
fprintf('Beamformer finished\n' )

trialsRPE = outcome_match.RPE(outcome_match.win~=0);
Rwin = trialsRPE>0; 
Rgamble = abs(trialsRPE)>5; % big gambles!
Rsafe = abs(trialsRPE)<3; % safe gambles

datapress= define_trials(buttonpress,data,bv_match.time,[-0.2 1]);
datapressc = cell2mat(datapress.trial);

ntrialsp = length(datapress.time);
nsamplesp = length(datapress.time{1});
VEpress(1:size(L,2)) = {0};
parfor ii = 1:length(L)
    if ~isempty(L{ii})
        w = W{ii};
        wnorm = w/sqrt( sum( (w*noise).^2) );

        dataloc = wnorm'*dataerfc;
        dataloc = reshape(dataloc, [nsamples, ntrials]);
       
        %         Pw{ii} = mean(dataloc(winlose==1,:),1);
        %         Pl{ii} = mean(dataloc(winlose==-1,:),1);
        VEp{ii} = mean(dataloc(:,Rwin),2);
        VEn{ii} = mean(dataloc(:,~Rwin),2);
%         VEg{ii} = mean(dataloc(:,Rgamble),2);
%         VEs{ii} = mean(dataloc(:,Rsafe),2);
        
        dataloc = wnorm'*datapressc;
        dataloc = reshape(dataloc, [nsamplesp, ntrialsp]);
        VEpress{ii} = mean(dataloc,2);
    end
    if mod(ii,300) == 0
        fprintf('SAM running %.1f perc.\n', ii/length(L)*100 )
    end
end

%%

tstep = 0.1;

for twind = -.2:tstep:0.6

VEerf = zeros(size(VEpress));
tt = dataerf.time{1}>twind & dataerf.time{1}<(twind+tstep) ; 
for ii = 1:length(L)
    if ~isempty(L{ii})
%         VEerf(ii) = mean(abs(VEpress{ii}(tt)) );
%         VEerf(ii) = max( abs(VEp{ii}(tt))  );
%         VEerf(ii) = max( (abs(VEp{ii}(tt))*nnz(Rwin) + abs(VEn{ii}(tt))*nnz(~Rwin))/length(Rwin)  );

        VE = (abs(VEp{ii})*nnz(Rwin) + abs(VEn{ii})*nnz(~Rwin))/length(Rwin);
        xx = (VE-mean(VE)) > 3*std(VE); 
        p = max(VE(tt & xx'));
        if isempty(p)
            VEerf(ii) = 0;
        else
            VEerf(ii) = p;
        end
        
    end
end

% Plot Beamformer
sourceout =[];
sourceout.pow = VEerf;
sourceout.dim = dataprep.leadfield.dim;
sourceout.inside = dataprep.leadfield.inside;
sourceout.pos  = dataprep.leadfield.pos;
cfg = [];
cfg.parameter = 'pow';
sourceout_Int  = ft_sourceinterpolate(cfg, sourceout , dataprep.mri);

crang = [0 10]*1e-3;
cfg = [];
cfg.method        = 'slice';
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap  = 'auto';
cfg.funcolorlim   = crang;
cfg.opacitylim = crang;
ft_sourceplot(cfg, sourceout_Int);

set(gcf,'name',sprintf('Win ERF, %.0f-%.0fms',twind*1000,(twind+tstep)*1000))
clc
end
%%
% x = [20 10 100
%     -60 20 80
%     -40 30 100
%     30 10 60
%     -50 -40 50];
    
x = [-10 -55 80
    40 -40 45
    40 40 45];  

x = [-45 -25 85
    -40 -30 50
    25  30 90
   ];

x  = [20 5  45];

voxpos = round(dataprep.leadfield.pos*1000); % convert to mm
figure
for ii = 1:size(x,1)
    
    xx = find(voxpos(:,1) == x(ii,1) & voxpos(:,2) == x(ii,2) & voxpos(:,3) == x(ii,3) );
    w = W{xx};
    dataloc = w'*dataerfc; % gamble outcome
    dataloc = reshape(dataloc, [nsamples, ntrials]); 
    
%     dataloc = w'*datapressc; % button presses
%     dataloc = reshape(dataloc, [nsamplesp, ntrialsp]); 
    VE = mean(dataloc,2)*1e9;
    if abs(max(VE)) > abs(min(VE))
        plot(dataerf.time{1}*1e3,VE)
    else
        plot(dataerf.time{1}*1e3,-VE)
    end
    
    hold on
end
% plot(dataerf.time{1},VEn{x})
legend(num2str(x))
xlabel('Time (ms)'); ylabel('nAm')
%% ERF analysis 
% Rsamples     = outcome_match.sample(outcome_match.win==1);
% dataerf = define_trials(Rsamples,data,bv_match.time,[-0.2 1]);
% 
% % dataerf = define_trials(buttonpress,data,bv_match.time,[-0.2 1]);
% avgpress = ft_timelockanalysis(cfg, dataerf);
% [coeff,score] = pca(avgpress.avg');
% 
% % make a topographic plot of the ERF, in steps of 5ms
% cfg = [];
% cfg.xlim = [0:6];
% cfg.colorbar = 'no';
% cfg.comment = '';
% cfg.showxlim = 'no';
% cfg.showzlim = 'no';
% cfg.layout = 'CTF275_helmet.mat';
% cfg.parameter = 'avg';
% 
% cfg.zlim = [ ] * 1e-13;
% figure
% avgpca = avgpress;
% avgpca.avg = coeff;
% avgpca.time = 1:size(coeff,2);
% avgpca.var = zeros(size(coeff));
% avgpca.dof = zeros(size(coeff));
% ft_topoplotER(cfg, avgpca);
% 
% figure
% for tt = 1:6
%    subplot(2,3,tt)
%    plot(avgpress.time,score(:,tt))
%    xlabel('time(s)'); ylabel('B(t), Tesla')
%    grid on
% end
%     
%     
% %%
% % make a topographic plot of the ERF, in steps of 50ms
% cfg = [];
% cfg.xlim = -0.1:0.05:1.15;
% cfg.colorbar = 'no';
% cfg.comment = '';
% cfg.showxlim = 'no';
% cfg.showzlim = 'no';
% cfg.layout = 'CTF275_helmet.mat';
% 
% cfg.zlim = [-1 1] * 1e-13;
% figure
% ft_topoplotER(cfg, avgoutc);
% 
% 
% %%
% Rsamples     = ITI_match.sample(ITI_match.sample~=0);
% dataerf = define_trials(Rsamples,data,bv_match.time,[-1 1]);
% 
% cfg = [];
% avgoutc = ft_timelockanalysis(cfg, dataerf);
% 
% 
% figure
% plot(1000*avgoutc.time, avgoutc.avg)  % convert time to ms
% xlabel('time (ms)')
% ylabel('field amplitude (T)')
% axis tight
% grid on
% 
% 
% figure
% plot(1000*avgpress.time, avgpress.avg')  % convert time to ms
% xlabel('time (ms)')
% ylabel('field amplitude (T)')
% axis tight
% grid on
% % xlim([-400 2000])
% %% It works with buttin presses!! Important to add symmetry
% cfg = [];
% cfg.latency = [0.2 0.4];  % specify latency window around M50 peak
% cfg.numdipoles = 2;
% cfg.symmetry = 'y';
% % cfg.hdmfile = 'SubjectBraille.hdm';
% cfg.mri = dataprep.mri;
% % cfg.grad            = sens;
% cfg.headmodel       = dataprep.headmodel;
% cfg.sourcemodel.unit   = 'cm';
% % cfg.feedback = 'textbar';
% cfg.resolution = 2;
% cfg.unit = 'cm';
% cfg.reducerank  =2;
% cfg.normalize  ='no';
% cfg.dipfit.checkinside  = true;
% % cfg.dip.pos     = [1.9 3.5 8.8];
% % cfg.gridsearch = 'no';
% 
% dipM50 = ft_dipolefitting(cfg, avgpress);
% 
% cfg = [];
% cfg.location = dipM50.dip.pos(1,:)*10;   % convert from cm to mm
% ft_sourceplot(cfg, dataprep.mri)