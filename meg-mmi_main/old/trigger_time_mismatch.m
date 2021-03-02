clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults

sample_mismatch = [];


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'];

for s = 3:4
    

sub = subn(s,:); % with pixel '24138'artefacts; '24103';    % no pixel '24172'; '24071'
data_path = ['/data/MBDU/MEG_MMI3/data/derivatives/sub-',sub,'/'];
cd(data_path)

data_name = [sub,'MMI_mmi3_proc.ds']; %Pre-processed at 1-150 Hz to adjust baseline

%% Read events



bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');
for ii = 1:length(bv_names)
    if strncmp(bv_names(ii).name,sub,5)
        bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
bv = readtable(bv_name,opts); % bahavioral data

ratemood_block = bv.blockHappyText_started;
ratemood_block(isnan(ratemood_block))= [] ;

event = ft_read_event(data_name);

% pix = ft_read_data(data_name,'chanindx','UADC016');

cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'no';
cfg.channel = 'UADC016';
cfg.demean = 'yes';
pix = ft_preprocessing(cfg);

pixm  = cell2mat(pix.trial)';
pix_white = pixm>2.5;
pix_sample = find(diff(pix_white)==1);

f =pix.fsample;

bv_answer = bv.fixCross_started;
bv_answer(isnan(bv_answer))=[];

bv_choice = bv.fixCross_2_started;
bv_choice(isnan(bv_choice)) = [];

bv_outcome = bv.fixCross_3_started;
bv_outcome(isnan(bv_outcome)) = [];

bv_mood_block = bv.blockHappyText_started;
bv_mood_block(isnan(bv_mood_block)) = [];

bv_slider_block = bv.blockHappySlider_started;
bv_slider_block(isnan(bv_slider_block)) = [];

bv_mood = bv.happyText_started;
bv_mood(isnan(bv_mood)) = [];

bv_slider = bv.happySlider_started;
bv_slider(isnan(bv_slider)) = [];

bv_rest = bv.endOfBlockText_started;
bv_rest(isnan(bv_rest)) = [];

% Check time difference between pixel and behavioral file
bv_all = sort([bv_answer; bv_choice; bv_outcome; bv_mood; bv_mood_block; bv_slider; bv_slider_block; bv_rest]);

bv_timeoffset = median(bv_all(1:length(pix_sample))-pix_sample/f); % difference in start time between the behavioral and MEG data.
if std(bv_all(1:length(pix_sample))-pix_sample/f) > 20e-3
    warning('MEG and behavioral file temporal mismatch > 20ms')
    figure; histogram(bv_all(1:length(pix_sample))-pix_sample/f);
    xlabel('Time off-set (seconds)');
end

answer_match = [];
choice_match = [];
outcome_match =[];
ITI_match = [];
mood_match = [];
blockmood_match = [];
slider_match = [];
blockslider_match = [];
% Find all choices
nn = 0;
for ii = 1:length(pix_sample)
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        answer_match.pix_index(nn) = ii;
        answer_match.bv_index(nn) = jj;
        answer_match.sample(nn) = pix_sample(ii);
        %        answer_match.RPE(nn) = bv.RPE(jj);
        %        answer_match.winAmount = bv.winAmount(jj);
    end
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_2_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        choice_match.pix_index(nn) = ii;
        choice_match.bv_index(nn) = jj;
        choice_match.sample(nn) = pix_sample(ii);
        if strcmp(bv.choice{jj},'gamble')
            choice_match.gamble(nn) = 1;
        elseif strcmp(bv.choice{jj},'certain')
            choice_match.gamble(nn) = 0;
        end
    end
    
    sdiff = abs(pix_sample(ii) - (bv.fixCross_3_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        outcome_match.pix_index(nn) = ii;
        outcome_match.bv_index(nn) = jj;
        outcome_match.sample(nn) = pix_sample(ii);
        if strcmp(bv.outcome{jj},'win')
            outcome_match.win(nn) = 1;
        elseif strcmp(bv.outcome{jj},'lose')
            outcome_match.win(nn) = -1;
        elseif strcmp(bv.outcome{jj},'certain')
            outcome_match.win(nn) = 0;
        end
    end
    
    % Take fixation cross at the end of each trial as baseline
    % lasts 2s at the end of each trial
    sdiff = abs(pix_sample(ii)-2*f - (bv.fixCross_ITI_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        ITI_match.pix_index(nn) = ii;
        ITI_match.bv_index(nn) = jj;
        ITI_match.sample(nn) = pix_sample(ii)-2*f;
        %             if strcmp(bv.outcome{jj},'win')
        %                 ITI_match.win(nn) = 1;
        %             elseif strcmp(bv.outcome{jj},'lose')
        %                 outcome_match.win(nn) = -1;
        %             elseif strcmp(bv.outcome{jj},'certain')
        %                 outcome_match.win(nn) = 0;
        %             end
    end
    
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.happySlider_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        slider_match.pix_index(nn) = ii;
        slider_match.bv_index(nn) = jj;
        slider_match.sample(nn) = pix_sample(ii);
        slider_match.mood(nn) = bv.happySlider_started(jj);
    end
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.blockHappySlider_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        blockslider_match.pix_index(nn) = ii;
        blockslider_match.bv_index(nn) = jj;
        blockslider_match.sample(nn) = pix_sample(ii);
        blockslider_match.mood(nn) = bv.blockHappySlider_started(jj);
    end
    
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.happyText_started-bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        mood_match.pix_index(nn) = ii;
        mood_match.bv_index(nn) = jj;
        mood_match.sample(nn) = pix_sample(ii);
        mood_match.mood(nn) = bv.happyText_started(jj);
    end
    
    % Rate mood
    sdiff = abs(pix_sample(ii) - (bv.blockHappyText_started -bv_timeoffset)*f);
    [mm, jj] = min(sdiff);
    if mm < .1*f
        nn = jj;
        blockmood_match.pix_index(nn) = ii;
        blockmood_match.bv_index(nn) = jj;
        blockmood_match.sample(nn) = pix_sample(ii);
        blockmood_match.mood(nn) = bv.blockHappyText_started(jj);
    end
    
    
end


win_sample = cell(size(event));
lose_sample = cell(size(event));
ratemood_sample = cell(size(event));
slider_sample = cell(size(event));
gamble_sample = cell(size(event));
certain_sample = cell(size(event));
answer_sample = cell(size(event));
rest_sample = cell(size(event));

for tt = 1:length(event)
    
    if strcmp(event(tt).type, 'UPPT001' ) % Slider appear 3s after        
        switch event(tt).value
            case 1
                answer_sample{tt} = event(tt).sample;
            case 2
                certain_sample{tt} = event(tt).sample;                
            case 4
                gamble_sample{tt} = event(tt).sample;
            case 8
                win_sample{tt} = event(tt).sample;
            case 16
                lose_sample{tt} = event(tt).sample;
            case 32
                ratemood_sample{tt} = event(tt).sample;
            case 64
                slider_sample{tt} = event(tt).sample;
            case 128
                rest_sample{tt} = event(tt).sample;
        end
    end
    
end

ratemood_sample = cell2mat(ratemood_sample);
lose_sample = cell2mat(lose_sample);
win_sample = cell2mat(win_sample);
gamble_sample = cell2mat(gamble_sample);
certain_sample = cell2mat(certain_sample);
answer_sample = cell2mat(answer_sample);
rest_sample = cell2mat(rest_sample);
slider_sample =cell2mat(slider_sample);


[~,~,m] = find(mood_match.sample);
[~,~,m1] = find(blockmood_match.sample);
sample_mismatch = cat(2,sample_mismatch ,sort([m,m1]) - ratemood_sample');

[~,~,m] = find(slider_match.sample);
[~,~,m1] = find(blockslider_match.sample);
sample_mismatch = cat(2,sample_mismatch ,sort([m,m1]) - slider_sample');

[~,~,m] = find(answer_match.sample);
sample_mismatch = cat(2,sample_mismatch ,m - answer_sample');


[~,~,m] = find(choice_match.sample);
if s == 3
m(31)= []; % missing trigger
end
sample_mismatch = cat(2,sample_mismatch ,m - sort([certain_sample;gamble_sample])');

m = choice_match.sample(logical(choice_match.gamble));
sample_mismatch = cat(2,sample_mismatch ,m - sort(gamble_sample)');

m = outcome_match.sample(logical(abs(outcome_match.win)));
sample_mismatch = cat(2,sample_mismatch ,m - sort([win_sample; lose_sample]'));

end

histogram(sample_mismatch)
% 24, points = 20ms is the mean/median delay
