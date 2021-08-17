% Analyzing LSD experiments

%%  we want to exclude fish that did not really sleep on night 1
% Not sure this makes sense, since we have mutants that could sleep very
% little; I will comment this out for now

% night1_sleep_slice = genox.sleep{1, 2}(841:1440, :); % get the sleep of night 1
% night1_sleep_sum = sum(night1_sleep_slice, 1); % sum that sleep for each fish

% good_sleepers = find(night1_sleep_sum > 5);
% number_of_good_sleepers = length(good_sleepers); % how many have we left?

% % find our baseline slice for SD part
% % Only using the first 4h here
% SD_baseline_slice = genox.sleep{1, 2}(841:1080, good_sleepers);
% % find our baseline slice for RS part
% RS_baseline_slice = genox.sleep{1, 2}(1201:1440, good_sleepers);
% % find the SD slice
% % Only using the first 4h here
% SD_slice = genox.sleep{1, 2}(2281:2520, good_sleepers);
% % find the RS slice
% RS_slice = genox.sleep{1, 2}(2641:2880, good_sleepers);

%%

% find our baseline slice for SD part; only using the first 4h here
SD_baseline_slice = genox.sleep{1, 2}(841:1080, :);
% find our baseline slice for RS part; only using the last 4h here
RS_baseline_slice = genox.sleep{1, 2}(1201:1440, :);
% find the SD slice; only using the first 4h here (total SD is 6hr)
SD_slice = genox.sleep{1, 2}(2281:2520, :);
% find the RS slice; last 4 hr
RS_slice = genox.sleep{1, 2}(2641:2880, :);


% for each fish, sum how much it slept during the SD_baseline_slice
SD_baseline_sum = sum(SD_baseline_slice, 1)
% for each fish, sum how much it slept during the RS_baseline_slice
RS_baseline_sum = sum(RS_baseline_slice, 1)
% for each fish, sum how much it slept during the SD phase
SD_sum = sum(SD_slice, 1)
% for each fish, sum how much it slept during the RS phase
RS_sum = sum(RS_slice, 1)

% we want to find the fish that lost at least some % of their sleep druing
% the SD phase. We compare sleep during SD with sleep during SD_baseline
% and if the index is above X, we keep those fish

SD_ratio = SD_sum./mean(SD_baseline_sum)
Sleep_lost_index = 1 - SD_ratio 
sleep_deprived_fish = find(Sleep_lost_index > 0.5)


% Re-filter for only SDF (Sleep Deprived Fish)
% find our baseline slice for SD part; only using the first 4h here
SD_baseline_slice_SDF = genox.sleep{1, 2}(841:1080, sleep_deprived_fish);
% find our baseline slice for RS part; only using the last 4h here
RS_baseline_slice_SDF = genox.sleep{1, 2}(1201:1440, sleep_deprived_fish);
% find the SD slice; only using the first 4h here (total SD is 6hr)
SD_slice_SDF = genox.sleep{1, 2}(2281:2520, sleep_deprived_fish);
% find the RS slice; last 4 hr
RS_slice_SDF = genox.sleep{1, 2}(2641:2880, sleep_deprived_fish);



% for each SD fish, sum how much it slept during the SD_baseline_slice
SD_baseline_sum_SDF = sum(SD_baseline_slice_SDF, 1)
% for each SD fish, sum how much it slept during the RS_baseline_slice
RS_baseline_sum_SDF = sum(RS_baseline_slice_SDF, 1)
% for each SD fish, sum how much it slept during the SD phase
SD_sum_SDF = sum(SD_slice_SDF, 1)
% for each SD fish, sum how much it slept during the RS phase
RS_sum_SDF = sum(RS_slice_SDF, 1)

% Having found the fish that were actually sleep deprived, we can use those
% to calculate what ratio of sleep was recovered after sleep deprivation.
% The general idea is:
% ratio_of_sleep_recovered = sleep_recovered/sleep_lost

% For the sleep_recovered I subtract the first 4 hr of light-out of the
% second day from the first day (correct?)
sleep_recovered = RS_sum_SDF - SD_baseline_sum_SDF
sleep_lost = SD_baseline_sum_SDF - SD_sum_SDF
ratio_of_sleep_recovered = (sleep_recovered./sleep_lost)'
recovery_sleep_norm = (RS_sum_SDF./SD_baseline_sum_SDF)'
recovery_sleep_norm_norm = (RS_sum_SDF/mean(SD_baseline_sum_SDF))'
bounce_effect = (RS_sum_SDF/mean(SD_sum_SDF))'

% Another idea is:
% SD_response = Recovery sleep (of each fish) / sleep_lost (either of each
% fish, or of the average for the whole genotype)

SD_response = (RS_sum_SDF./sleep_lost)'
SD_response_norm = (RS_sum_SDF./mean(sleep_lost))' % this is the relevant one






%% Random stuff

% Mean_Sleep_lost_index = 1 - mean(SD_ratio)
% 
% 
% RS_baseline_mean = mean(RS_baseline_sum)
% 
% RS_difference = RS_sum - RS_baseline_mean
% 
% Sleep_lost = sum(SD_baseline_sum) - sum(SD_sum)
% Mean_Sleep_recovered_1 = sum(RS_sum) - sum(RS_baseline_sum)
% Mean_Sleep_recovered_2 = sum(RS_sum) - sum(SD_baseline_sum) % this seems to be the better way to do this
% 
% 
% Mean_Percent_recovered_1 = 100*Mean_Sleep_recovered_1/Sleep_lost
% 
% Mean_Percent_recovered_2 = 100*Mean_Sleep_recovered_2/Sleep_lost
% 
% Mean_Sleep_lost_index


% 
% 
% RS_baseline_mean_SDF = mean(RS_baseline_sum_SDF)
% 
% RS_difference_SDF = RS_sum_SDF - RS_baseline_mean_SDF
% 
% Sleep_lost_SDF = SD_baseline_sum_SDF - SD_sum_SDF
% % Mean_Sleep_recovered_1_SDF = sum(RS_sum_SDF) - sum(RS_baseline_sum_SDF)
%  Mean_Sleep_recovered_2_SDF = sum(RS_sum_SDF) - sum(SD_baseline_sum_SDF) % this seems to be the better way to do this
% 
% Sleep_recovered_2_SDF = RS_sum_SDF - SD_baseline_sum_SDF
% 
% ratio_recovered_2_SDF = Sleep_recovered_2_SDF./Sleep_lost_SDF
% 
% % Mean_Percent_recovered_1_SDF = 100*Mean_Sleep_recovered_1_SDF/Sleep_lost_SDF
% % Mean_Percent_recovered_2_SDF = 100*Mean_Sleep_recovered_2_SDF/Sleep_lost_SDF

