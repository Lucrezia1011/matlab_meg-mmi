function [R,Rlta]= mmiBehavioralAll(sub)

bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');
for ii = 1:length(bv_names)
    if strncmp(bv_names(ii).name,sub,5)
        bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
bv = readtable(bv_name,opts); % bahavioral data

indGamble = strcmp(bv.outcome,'win') | strcmp(bv.outcome,'lose');

EltaH = cumsum(bv.outcomeAmount(~isnan(bv.outcomeAmount)))./(1:81)'; % Expectation, defined by Hanna
RltaH = bv.outcomeAmount(~isnan(bv.outcomeAmount)) - EltaH;
Elta = zeros(size(EltaH));
Elta(2:end) = EltaH(1:end-1);
Rlta = bv.outcomeAmount(~isnan(bv.outcomeAmount)) - Elta;
RPE = bv.RPE(~isnan(bv.outcomeAmount));
indGamble = indGamble(~isnan(bv.RPE));

R = RPE(indGamble)/std(RPE(indGamble));
Rlta = Rlta(indGamble)/std(Rlta(indGamble));