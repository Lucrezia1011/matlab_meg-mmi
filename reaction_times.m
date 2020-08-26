
meginfo = readtable('~/MEG_participantsinfo.csv');

RT = struct;
RT.reactionTime = [];
RT.choice = [];
RT.block = [];

Nlistx = [1:9,11:56];  % only exclude subject 10
for sn = Nlistx
    sub = num2str(meginfo.SDAN(sn));


    bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');
    for ii = 1:length(bv_names)
        if strncmp(bv_names(ii).name,sub,5)
            bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
        end
    end

    opts = detectImportOptions(bv_name);
    bv = readtable(bv_name,opts); % bahavioral data
    
    RT.reactionTime = cat(1,RT.reactionTime,bv.choiceKey_rt(12+(1:81))); % reacion time
    RT.choice = cat(1,RT.choice,bv.choice(12+(1:81)));
    RT.block = cat(1,RT.block,bv.blockNum(12+(1:81)));

end

figure; set(gcf,'color','w')
subplot(131); histogram(RT.reactionTime)
xlabel('reaction time (s)')
subplot(132); boxplot(RT.reactionTime,RT.choice)
xlabel('choice'); ylabel('reaction time (s)')
subplot(133); boxplot(RT.reactionTime,RT.block)
xlabel('block'); ylabel('reaction time (s)')

[h,p] = kstest2(RT.reactionTime(strcmp(RT.choice,'gamble')),RT.reactionTime(strcmp(RT.choice,'certain')));

RTblocks = [RT.reactionTime(RT.block==1),RT.reactionTime(RT.block==2),RT.reactionTime(RT.block==3)];
% [iia,iib]=find(isnan(RTblocks));
% RTblocks(iia,:) = [];
[p,tbl,stats] = anova1(RTblocks);
% 

