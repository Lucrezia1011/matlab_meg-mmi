
meginfo = readtable('~/MEG_participantsinfo.csv');

RT = struct;
RT.reactionTime = [];
RT.choice = [];
RT.block = [];
mood = NaN(81,56);

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

    mood(:,sn) = bv.happySlider_response(12+(1:81));
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
tt =~isnan(mood(:,1));
trials = 1:81;
figure; clf;set(gcf,'color','w')
plot(trials(tt),mean(mood(tt,Nlistx),2),'k','LineWidth',3)
hold on
fill(27+[0 27 27 0],[0 0 1 1],[0 0 0.5],'facealpha',0.1,'edgecolor','none')
fill([trials(tt) fliplr(trials(tt))],...
    [mean(mood(tt,Nlistx),2)+std(mood(tt,Nlistx),0,2)/sqrt(length(Nlistx)); ...
    flipud(mean(mood(tt,Nlistx),2)-std(mood(tt,Nlistx),0,2)/sqrt(length(Nlistx)))],...
    [0.5 0.5 0.5],'facealpha',0.4,'edgecolor','none')
xlabel('Trial'); ylabel('Mood')
title(sprintf('Average mood over all participants, N=%.0f',nnz(Nlistx)))
xlim([3 81]); ylim([0.45 0.8])


hvs = intersect(Nlistx,find(strcmp(meginfo.Group,'HV')));
mdds = intersect(Nlistx,find(strcmp(meginfo.Group,'MDD')));
figure; clf;set(gcf,'color','w')

subplot(121)
plot(trials(tt),mean(mood(tt,hvs),2),'k','LineWidth',3)
hold on
fill(27+[0 27 27 0],[0 0 1 1],[0 0 0.5],'facealpha',0.1,'edgecolor','none')
fill([trials(tt) fliplr(trials(tt))],...
    [mean(mood(tt,hvs),2)+std(mood(tt,hvs),0,2)/sqrt(length(hvs)); ...
    flipud(mean(mood(tt,hvs),2)-std(mood(tt,hvs),0,2)/sqrt(length(hvs)))],...
    [0.5 0.5 0.5],'facealpha',0.4,'edgecolor','none')
xlabel('Trial'); ylabel('Mood')
title(sprintf('Average mood over all HVs, N=%.0f',nnz(hvs)))
xlim([3 81]); ylim([0.4 0.9])

subplot(122)
plot(trials(tt),mean(mood(tt,mdds),2),'k','LineWidth',3)
hold on
fill(27+[0 27 27 0],[0 0 1 1],[0 0 0.5],'facealpha',0.1,'edgecolor','none')
fill([trials(tt) fliplr(trials(tt))],...
    [mean(mood(tt,mdds),2)+std(mood(tt,mdds),0,2)/sqrt(length(mdds)); ...
    flipud(mean(mood(tt,mdds),2)-std(mood(tt,mdds),0,2)/sqrt(length(mdds)))],...
    [0.5 0.5 0.5],'facealpha',0.4,'edgecolor','none')
xlabel('Trial'); ylabel('Mood')
title(sprintf('Average mood over all MDDs, N=%.0f',nnz(mdds)))
xlim([3 81]); ylim([0.4 0.9])


figure; clf;set(gcf,'color','w')
for n = 1:18
    subplot(4,5,n)
    plot(trials(tt),mood(tt,Nlistx((n-1)*3+1:n*3)))
    ylim([0 1])
    title(sprintf('subjects %.0f,%.0f,%.0f,%.0f',Nlistx((n-1)*3+1:n*3)))
end
n=19;
subplot(4,5,n)
plot(trials(tt),mood(tt,Nlistx((n-1)*3+1:end)))
ylim([0 1])
title(sprintf('subjects %.0f,%.0f,%.0f,%.0f',Nlistx((n-1)*3+1:end)))