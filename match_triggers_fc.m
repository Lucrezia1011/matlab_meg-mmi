function [bv_match,bv] = match_triggers_fc(data_name)
% f = resampling frequency

sub = data_name(1:5);

%% Read behavioral file
bv_names = dir('/data/MBDU/MEG_MMI3/data/behavioral/');
for ii = 1:length(bv_names)
    if strncmp(bv_names(ii).name,sub,5)
        bv_name = ['/data/MBDU/MEG_MMI3/data/behavioral/',bv_names(ii).name];
    end
end

opts = detectImportOptions(bv_name);
bv = readtable(bv_name,opts); % bahavioral data


%% Select data


restt_found = 0;

disp(data_name)

% Read events: Double check for non proc data

% pix = ft_read_data(data_name,'chanindx','UADC016');
% Outdated way to read light pixel channel
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = 'SCLK01'; % Time channel!
time = ft_preprocessing(cfg);
f =time.fsample;

% read LIGHT marker
ii = 0;
event = ft_read_event(data_name);
pix_sample = zeros(size(event));
for tt = 1:length(event)
    if strcmp(event(tt).type, 'LIGHT' )
        ii = ii +1;
        pix_sample(ii) = event(tt).sample;
    end
end
d = diff(pix_sample);
l = 0.25*f; % Eliminates triggers separated by less than 250ms (e.g. double marking of single trigger)
ii = find(d<l);
pix_sample(ii+1) = []; % correctly identifies triggers and eliminates zeros


trig_sample = zeros(size(event,1),2);
ii = 0;
for tt = 1:length(event)
    if strcmp(event(tt).type, 'UPPT001' )
        ii = ii +1;
        trig_sample(ii,1) = event(tt).sample + 24; % Average 24 sample delay
        trig_sample(ii,2) = event(tt).value;
    end
end
trig_sample = trig_sample(1:ii,:);

% Uses trigger samples when pixel is missing
if nnz(pix_sample) == 0
    warning('No LIGHT trigger detected: switching to UPPT001')
    pix_sample = trig_sample(:,1);
else % Checks if any light triggers are missing (e.g. when participants respons faster than 200ms)
    % for participant 24199 (pressing gamble option too quickly) trig
    % value is = 5 instead of normal 4
    d = diff(trig_sample);
    ii = find(d(:,1)<l & d(:,2)~= 0); % finds different triggers closer than 250ms
    
    for jj = ii'
        sdiff = abs(trig_sample(jj,1) - pix_sample);
        [mm,kk] = min(sdiff);
        pix_sample(end+1) = trig_sample(jj+1,1)+mm;
    end
    pix_sample = sort(pix_sample);
    
end

% Actual recording time (account for skipped trials)
pix_time = time.trial{1}(pix_sample); %  - time.trial{1}(1)
% 
% bv_answer = bv.fixCross_started;
% bv_answer(isnan(bv_answer))=[];
% 
% bv_choice = bv.fixCross_2_started;
% bv_choice(isnan(bv_choice)) = [];
% 
% bv_outcome = bv.fixCross_3_started;
% bv_outcome(isnan(bv_outcome)) = [];
% 
bv_mood_block = bv.blockHappyText_started;
bv_mood_block(isnan(bv_mood_block)) = [];
% 
% bv_slider_block = bv.blockHappySlider_started;
% bv_slider_block(isnan(bv_slider_block)) = [];
% 
% bv_mood = bv.happyText_started;
% bv_mood(isnan(bv_mood)) = [];
% 
% bv_slider = bv.happySlider_started;
% bv_slider(isnan(bv_slider)) = [];

bv_rest = bv.endOfBlockText_started;
bv_rest(isnan(bv_rest)) = [];
%%

% Check time difference between pixel and behavioral file
%     bv_all = sort([bv_answer; bv_choice; bv_outcome; bv_mood; bv_mood_block; bv_slider; bv_slider_block; bv_rest]);
%     if length(bv_all) >= length(pix_sample) % Checks there are no extra triggers due to slider pixel
%     bv_timeoffset = median(bv_all(1:length(pix_sample))-pix_time); % difference in start time between the behavioral and MEG data.
%         if std(bv_all(1:length(pix_sample))-pix_sample/f) > 20e-3
%             warning('MEG and behavioral file temporal mismatch > 20ms')
%             figure; histogram(bv_all(1:length(pix_sample))-pix_time);
%             xlabel('Time off-set (seconds)');

%             warning('Using first trigger as time offset')
%             bv_timeoffset = bv_all(1) - pix_time(1);

% Aligns behavioral file and MEG with through rest trigger
rest_sample = [];
for tt = 1:length(event)
    if strcmp(event(tt).type, 'UPPT001' ) && event(tt).value == 128
        rest_sample = [rest_sample, event(tt).sample + 24]; % Average 24 sample delay
    end
end


if ~isempty(rest_sample)
    % Checks all combinations of rest triggers in MEG and
    % behavioral files.
    rest_diff = bv_rest - time.trial{1}(rest_sample);
    if length(rest_sample) == 3
        restt_found = restt_found + 3;
        disp('Matching all 3 rest trials to behavioral file')
        bv_timeoffset = mean(diag(rest_diff));
    elseif length(rest_sample) == 2
        restt_found = restt_found+2;
        d = diag(rest_diff);
        if abs(diff(d)) > 0.1 %
            d = diag(rest_diff,-1);
            disp('Matching 2^ and 3^ rest triggers to behavioral file')
        else
            disp('Matching 1^ and 2^ rest triggers to behavioral file')
        end
        bv_timeoffset  = mean(d);
    elseif length(rest_sample) == 1
        restt_found = restt_found +1;
        disp(['Found only one rest trigger: assuming to be ',num2str(restt_found),'^'])
        bv_timeoffset  = rest_diff(restt_found);
    end
else
    if strcmp(data_name(19),'3')
        bv_timeoffset = bv_mood_block(3) - pix_time(1);
        warning('No rest triggers found: matching first trigger with third block mood rating')
    else
        warning('No rest triggers found: matching first trigger with first block mood rating')
        bv_timeoffset = bv_mood_block(1) - pix_time(1);
    end
end

match_check = 1;

try
    while match_check
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
        
        w = .1;
        for ii = 1:length(pix_sample)
            
            sdiff = abs(pix_time(ii) - (bv.fixCross_started - bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                answer_match.pix_index(nn) = ii;
                answer_match.bv_index(nn) = jj;
                answer_match.sample(nn) = pix_sample(ii);
                answer_match.choice(nn) = bv.choiceKey_keys(jj);
                answer_match.RT(nn)  =  bv.choiceKey_rt(jj);
                %        answer_match.RPE(nn) = bv.RPE(jj);
                %        answer_match.winAmount = bv.winAmount(jj);
            end
            
            sdiff = abs(pix_time(ii) - (bv.fixCross_2_started- bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                choice_match.pix_index(nn) = ii;
                choice_match.bv_index(nn) = jj;
                choice_match.sample(nn) = pix_sample(ii);
                if strcmp(bv.choice{jj},'gamble')
                    choice_match.gamble(nn) = 1;
                elseif strcmp(bv.choice{jj},'certain')
                    choice_match.gamble(nn) = 0;
                end
                choice_match.choice(nn) = bv.choiceKey_keys(jj);
            end
            
            sdiff = abs(pix_time(ii) - (bv.fixCross_3_started-bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
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
                outcome_match.winamount(nn) = bv.winAmount(jj);
                outcome_match.loseamount(nn) = bv.loseAmount(jj);
                outcome_match.RPE(nn) = bv.RPE(jj);

                
            end
            
            % Take fixation cross at the end of each trial as baseline
            % lasts 2s at the end of each trial
            sdiff = abs(pix_time(ii)-2 - (bv.fixCross_ITI_started-bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
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
            sdiff = abs(pix_time(ii) - (bv.happySlider_started-bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                slider_match.pix_index(nn) = ii;
                slider_match.bv_index(nn) = jj;
                slider_match.sample(nn) = pix_sample(ii);
                slider_match.mood(nn) = bv.happySlider_response(jj);
            end
            
            % Rate mood
            sdiff = abs(pix_time(ii) - (bv.blockHappySlider_started -bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                blockslider_match.pix_index(nn) = ii;
                blockslider_match.bv_index(nn) = jj;
                blockslider_match.sample(nn) = pix_sample(ii);
                blockslider_match.mood(nn) = bv.blockHappySlider_response(jj);
            end
            
            
            % Rate mood
            sdiff = abs(pix_time(ii) - (bv.happyText_started-bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                mood_match.pix_index(nn) = ii;
                mood_match.bv_index(nn) = jj;
                mood_match.sample(nn) = pix_sample(ii);
                mood_match.mood(nn) = bv.happySlider_response(jj);
            end
            
            % Rate mood
            sdiff = abs(pix_time(ii) - (bv.blockHappyText_started -bv_timeoffset));
            [mm, jj] = min(sdiff);
            if mm < w
                nn = jj;
                blockmood_match.pix_index(nn) = ii;
                blockmood_match.bv_index(nn) = jj;
                blockmood_match.sample(nn) = pix_sample(ii);
                blockmood_match.mood(nn) = bv.blockHappySlider_response(jj);
            end
            
            
        end
        
        % Did not match triggers with behavioral file
        if isempty(answer_match) || isempty(outcome_match) || ...
                isempty(choice_match) || isempty(mood_match) || isempty(slider_match) || nnz(mood_match.sample)<5
            
            if ~isempty(rest_sample)
                % If rest trigger is present, shift by one
                if length(rest_sample) == 1 && restt_found < 3
                    restt_found = restt_found +1;
                    disp(['Found only one rest trigger: assuming to be ',num2str(restt_found),'^'])
                    bv_timeoffset  = rest_diff(restt_found);
                else % if 2 or 3 rest triggers are present matching should not be ambivalent
                    error('Triggers not matching behavioral file!')
                end
            else
                error('No rest triggers found: Triggers not matching behavioral file!')
                %             disp('No rest triggers found: matching first trigger with first block mood rating')
                %             bv_timeoffset = bv_mood_block(1) - pix_time(1);
            end
            
        else
            match_check = 0;
        end
        
    end
end

%% Button presses
cfg = [];
cfg.dataset = data_name;
cfg.continuous = 'yes';
cfg.channel = {'UADC005';'UADC006';'UADC007'}; % Time channel!
buttons = ft_preprocessing(cfg);

buttonsd = diff(buttons.trial{1}');
buttonpress = buttonsd>1.5;
[samples,~] = find(buttonpress);
samples = sort(samples);

samples(find(diff(samples)<0.2*f)+1) = [];
%%

bv_match = struct;
% bv_match.deltrial = find(diff(time.trial{1})>1);
bv_match.time = time.trial{1};
bv_match.answer = answer_match;
bv_match.choice = choice_match;
bv_match.outcome = outcome_match;
bv_match.ratemood = mood_match;
bv_match.blockmood = blockmood_match;
bv_match.slider = slider_match;
bv_match.blockslider = blockslider_match;
bv_match.ITI = ITI_match;
bv_match.buttonpress = samples;
end

