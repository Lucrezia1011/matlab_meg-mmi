function LTAvars = LTA_calc(bv)


% Standard model
A = bv.outcomeAmount; % Outcome
A(isnan(A)) = [];
ntrials = length(A);

% Only for gamble, to check for outome window, anticioation and
% consumation
Es = (bv.winAmount + bv.loseAmount )/2;
Es(isnan(Es)) = [];
Rs = A - Es;

% LTA model
Elta = cumsum(A)./(1:ntrials)'; % Expectation, defined by Hanna
Elta(2:end) = Elta(1:end-1); 
Elta(1) = 0; % or Es(1)?
Rlta = A - Elta; % Assume RPE of first trial is A


bestfit_name = '/data/MBDU/MEG_MMI3/data/behavioral/closed_LTA_coefs.csv';
opts = detectImportOptions(bestfit_name);
bf_pars = readtable(bestfit_name,opts); 

bestfit_sub = bf_pars(bf_pars.Var1 == bv.participant(1),:);

g = bestfit_sub.gamma;
% g = .8;

sumE = zeros(ntrials,1);
sumR = zeros(ntrials,1);

for t = 1:ntrials
    sumE(t) = sum( g.^(0:(t-1))' .* Elta(t:-1:1) );
    sumR(t) = sum( g.^(0:(t-1))' .* Rlta(t:-1:1) );
end


bestfit_mood = bestfit_sub.m_0 + bestfit_sub.err_r + ...
    bestfit_sub.w_LTE*sumE + bestfit_sub.w_LTR*sumR;

LTAvars = struct;
LTAvars.E = Es;
LTAvars.R =  Rs;
LTAvars.E_LTA  = Elta;
LTAvars.E_sum = sumE;
LTAvars.R_LTA  = Rlta;
LTAvars.R_sum = sumR;
LTAvars.M = bestfit_mood;
LTAvars.M0 = bestfit_sub.m_0;
LTAvars.betaE = bestfit_sub.w_LTE;
LTAvars.betaR = bestfit_sub.w_LTR;

% tmood= bv.trialNumber(~isnan(bv.happySlider_response))+1;
% mood = bv.happySlider_response(tmood);
% figure
% plot(tmood,mood,'-o')
% hold on
% bestfit_mood = bestfit_sub.m_0 + bestfit_sub.err_r + ...
%     bestfit_sub.w_LTE*sumE(tmood) + bestfit_sub.w_LTR*sumR(tmood);
% plot(tmood,bestfit_mood,'-o')

