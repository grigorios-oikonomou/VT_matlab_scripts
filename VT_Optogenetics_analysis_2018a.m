%{
NOTES
> Script for the analysis of optogenetic experiments integrated at 1min intervals
> Grigorios Oikonomou - first version 2017
> 2018a version also graphs the 30 min after the end of the trial
I have marked the changes with  *** 2018a change ***
%}

%% Housekeeping
tic
close all, clear all; % clean up
working_path = '/Users/grigoriosoikonomou/Documents/Prober_Lab/Videotrackers/analysis_software/analysis2/';
%% Define the number, duration and timing of the Trials

trials = [822,912,1002,1092,1182,1272]; % night trials
% *** 2018a change *** to keep the trials store in the genox
genox.trials = trials; 

number_of_trials = length(trials);
length_of_trials = 30; % in minutes since we are doing minute analysis
for k=1:number_of_trials
    PSpeakT{k} = trials(k);
end

%{
PS is for Post Stimulus
BL is for BaseLine
GR is for GRaphing
PSpeakTx is the timepoint (minute) with the highest peak
The 3 min after the peak (3 min after LED-on) are excluded from the 
analysis (but are visualized)
The 2 min before Peak+Trial duration are also excluded from analysis, but
visualized (approximation of LED-off)
%}

for k=1:number_of_trials
    PSstartT{k} = PSpeakT{k} + 3; % the start of the PS is 3 minutes after the peak of the light-on response
    PSendT{k} = PSpeakT{k} + length_of_trials - 2; % the end of the PS is 2 minutes before the end trial (light-off response)
    BLendT{k} = PSpeakT{k} - 3; % the end of the BL is 3 minutes before the peak of the light-on response
    BLstartT{k} = PSpeakT{k} - length_of_trials + 2; % the start of the BL is set so that BL = PS (converse of PSendT)
    % *** 2018a change ***
    %graph_time_end{k} = PSendT{k}; % graph_time_end{k} = PSpeakT{k} + length_of_trials;
    %graph_time_start{k} = BLstartT{k}; % graph_time_start{k} = PSpeakT{k} - length_of_trials;
    graph_time_start{k} = PSpeakT{k} - length_of_trials; % graph_time_start{k} = PSpeakT{k} - length_of_trials;
    graph_time_end{k} = PSpeakT{k} + 2*length_of_trials; % graph_time_end{k} = PSpeakT{k} + length_of_trials;
end
%% Import genotype list(listname) and dataset

listname = '180709_4B_genotype_2.txt'; % the name of the genotype file
filename = '180709_4B_DATA.txt'; % the name of the data file

gridmap = importdata(strcat(working_path,'_genotype_lists/', listname), '\t', 2);
dataset = importdata(strcat(working_path,'_data_files/', filename), '\t', 2);

%% Import the names of the genotypes and sort the data in different cells based on genotypes 

% determine the number of genotypes you have
number_of_genos = length(gridmap.data(1,:));


% Split the data by genotypes. The genox.data{i} grabs for each genotype i
% the data from the full set dataset.data
for i=1:number_of_genos
    genox.name{i} = gridmap.colheaders{i};
    genox.data{i} = dataset.data(:, gridmap.data(~isnan(gridmap.data(:,i)),i)+2);
end


%% Save the fish numbers
% make a cell with doubles that are columns of the genotype numbers 
for i=1:number_of_genos
    genox.fishID{i} = gridmap.data(~isnan(gridmap.data(:,i)), i);
end


%% Create 60sec activity data (in case of 10sec integration)
% for i=1:number_of_genos
%     for j=6:6:length(genox.activity{i})
%         genox.activity60sec{i}(j/6,:) = sum(genox.data{i}(j-5:j, :))
%     end
% end

%% For 60 sec integrated data, just copy the data
for i=1:number_of_genos
            genox.activity60sec{i} = genox.data{i};
    
end

%% Clean up noise
% % You do not need this; you leave it in for activity, but for sleep, it is
% % cleaned up at a later step; I will keep this here in case the data gets
% % messy at some point in the future
% for i = 1:number_of_genos
%     genox.activity60sec{i}(genox.data{i} <= 0.1 ) = 0; %0.1 is noise for 60 sec integration; if equal or below that make it 0
% end

%% Get the genox.time60sec
% For 1min integrated data, just count the rows of the dataset 
for j=1:length(dataset.data)
    genox.time60sec(j) = (j);
end

% % For 10sec integrated data:
% % First count the 10sec intervals
% for j=1:length(genox.lightschedule)
%     genox.time10sec(j) = (j)
% end
% % count the 60sec (1 minute) intervals (1/1 of the total since we did 1min
% % integration)
% for j=1:floor(length(genox.lightschedule)/1)
%     genox.time60sec(j) = (j)
% end

%% ***************************   Stuff from original Jason script that I do not need  ***************************

% % Get the zeitgeber info from the dataset file (last column)
% % Just copy all the rows of the last column
% 
% genox.zeitgeber = dataset.data(:, end)
% 
% % Get the light schedule info from the zeitgeber info
% % lights are on at time 0, lights are off at time 14
% 
% genox.lightschedule = genox.zeitgeber % first just copy the zeitgeber
% genox.lightschedule(find(genox.lightschedule < 14)) = 0 % if <14 then lights are on, therefore 0
% genox.lightschedule(find(genox.lightschedule >= 14)) =1 % if >= 14 lights are off, therefore 1
% 
% % Remember: the stamp is 0 for day, and 1 for night.

%% Get Mean activity for each geno and graph that

%data
for i = 1:number_of_genos %loop through all genos
   % genox.activityMean60sec{i} = (sum(2, genox.activitySum10sec{i}(j-5:j))/length(genox.fishID{i}))' % for 10 sec integration
   % genox.activityMean60sec{i} = sum(genox.activity60sec{i}, 2)/length(genox.fishID{i}); % or just = mean(genox.activity60sec{i}, 2)
   genox.activityMean60sec{i} = mean(genox.activity60sec{i}, 2);
end

%graph
figure 
hold on
xlabel('Time since start (min)')
ylabel('Activity (sec/min)')
title (strcat('Mean Fish Activity (1min bins) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

for i=1:number_of_genos
    plot (genox.time60sec, genox.activityMean60sec{i})
    legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% axis and legend
legend(legendname)
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), '_Opto_Activity.fig'))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), '_Opto_Activity.tif'])

%% Create 60sec sleep data

for i=1:number_of_genos
    genox.sleep{i}=genox.activity60sec{i}; % first just copy the activity data
    genox.sleep{i}(genox.activity60sec{i} <= 0.1 ) = 1; % 0.1 is the minimum value given by the VT; if equal or below that, they have not moved, so they are asleep
    genox.sleep{i}(genox.activity60sec{i} > 0.1 ) = 0; % if above that, movement, hence 0 sleep  
end

%% Calculate and graph mean sleep

for i=1:number_of_genos
    genox.sleepMean{i} = sum(genox.sleep{i}, 2)/length(genox.fishID{i}); % or just = mean(genox.sleep{i}, 2)
    genox.sleepMean{i} = mean(genox.sleep{i}, 2);
end

%graph
figure
hold on
xlabel('Time since start (min)')
ylabel('Ratio of fish sleeping')
title (strcat('Fish sleep ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

for i=1:number_of_genos
    plot (genox.time60sec, genox.sleepMean{i})
    legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% axis and legend
legend(legendname)
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), '_Opto_Sleep.fig'))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), '_Opto_Sleep.tif'])

%% **************************************************    Activity and Normalized activity    **************************************************    
genox.SideBySideActivityNormal_All_T = [];

for k=1:number_of_trials
    
    % First get the BL and PS and GR activity regions
    for i = 1:number_of_genos %loop through all genos
        genox.BLactivityT{k}{i} = genox.activity60sec{i}(BLstartT{k}:BLendT{k},:); % BLactivity for all fist of this genotype in this trial; columns are fish, rows are minutes
        genox.PSactivityT{k}{i} = genox.activity60sec{i}(PSstartT{k}:PSendT{k},:); % PSactivity for all fist of this genotype in this trial; columns are fish, rows are minutes
        % *** 2018a change ***
        %genox.GRactivityT{k}{i} = genox.activity60sec{i}(GRstartT{k}:GRendT{k},:); % PSactivity for all fist of this genotype in this trial; columns are fish, rows are minutes
    end
    
    % Mean of Mean BL activity for each geno 
    for i = 1:number_of_genos %loop through all genos
        % genox.BLactivityMeanOfMeanT{k}{i} = mean((sum(genox.BLactivityT{k}{i})/(BLendT{k} - BLstartT{k} +1)), 2); % Mean BL activity for each fish, then mean that (mean of means)
        genox.BLactivityMeanOfMeanT{k}{i} = mean2(genox.BLactivityT{k}{i}); % mean2 is awesome!
    end
    
    % Mean PS activity for each fish; Means for all fish
    for i = 1:number_of_genos %loop through all genos
        % genox.PSactivityMeanT{k}{i} = sum(genox.PSactivityT{k}{i})/(PSendT{k} - PSstartT{k} +1); % mean PS activity for each fish
        genox.PSactivityMeanT{k}{i} = mean(genox.PSactivityT{k}{i}); % mean PS activity for each fish
    end
    
    % Divide the Mean PS activity of each fish with the mean BL for the genotype
    % to get the normalized PS activity
    for i = 1:number_of_genos %loop through all genos
        genox.PS_BLNormalT{k}{i} = (genox.PSactivityMeanT{k}{i}/genox.BLactivityMeanOfMeanT{k}{i})';
    end
           
    genox.SideBySideActivityNormalT{k} = padcat(genox.PS_BLNormalT{k}{1}, genox.PS_BLNormalT{k}{2});
    genox.SideBySideActivityNormal_All_T = cat(1, genox.SideBySideActivityNormal_All_T, genox.SideBySideActivityNormalT{k});
    
end


%% Activity graph

% Get all the trial activity slices next to each other
genox.all_trials_activity{i} = [];
for i = 1:number_of_genos
    for k = 1:number_of_trials
        previous_slices = genox.all_trials_activity{i};
        current_slice = genox.activity60sec{i}(graph_time_start{k}:graph_time_end{k}, :);
        genox.all_trials_activity{i} = cat(2, previous_slices, current_slice);
    end
end

% Normalize with graph_BL activity
genox.all_trials_activity_normalized{i} = [];
genox.all_trials_mean_activity_normalized{i} = [];
for i = 1:number_of_genos
    genox.graph_BL_activity{i} = mean2(genox.all_trials_activity{i}(3:(length_of_trials - 2), :));
    genox.all_trials_activity_normalized{i} = genox.all_trials_activity{i}/genox.graph_BL_activity{i};
    genox.all_trials_mean_activity_normalized{i} = mean(genox.all_trials_activity{i}/genox.graph_BL_activity{i}, 2);
end

% Make figures

% Activity figure with no error bars
figure
hold on

xlabel('Time (min)')
ylabel(strcat('Normalized Average Activity (a.u.)'))
title (strcat('Normalized Average Fish Activity (All Trials) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

% *** 2018a change ***
% x_axis = ((-length_of_trials + 2):(length_of_trials-2)); % exclude 2min from either side
x_axis = ((-length_of_trials):(2 * length_of_trials)); % 

for i=1:number_of_genos
    plot(x_axis, genox.all_trials_mean_activity_normalized{i})
    legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% axis and legend
legend(legendname, 'Location','NorthWest')
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Average_Activity_All_Trials', '.fig')))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials', '.tif')])
% save svg
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dsvg', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials', '.svg')])


% Activity Figure with shaded error bars
figure
hold on

xlabel('Time (min)')
ylabel(strcat('Normalized Average Activity (a.u.)'))
title (strcat('Normalized Average Fish Activity (All Trials) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

% x_axis = ((-length_of_trials + 2):(length_of_trials-2)); % exclude 2min from either side
x_axis = ((-length_of_trials):(2*length_of_trials)); % *** 2018a change ***

for i=1:number_of_genos
    
    if i == 1
        c = [0 0 1]; %blue
        s = 'b';
    elseif (i == 2 && number_of_genos == 2)
        c =  [1 0 0]; %red
        s = 'r';
    elseif (i == 2 && number_of_genos > 2)
        c = [0 1 1]; %cyan
        s = 'c';
    elseif (i == 3)
        c =  [1 0 0]; %red
        s = 'r';
    elseif i == 4
        c = [1 0 1]; %magenta/purple
        s = 'm';
    end
        
    sem = (nanstd(genox.all_trials_activity_normalized{i}, 0, 2)/sqrt(size(genox.all_trials_activity_normalized{i},2)))';
    shadedErrorBar(x_axis, genox.all_trials_mean_activity_normalized{i}, sem, 'lineprops', s);
    % legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% % axis and legend
% legend(legendname, 'Location','NorthWest')
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Average_Activity_All_Trials_shaded', '.fig')))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials_shaded', '.tif')])
% save svg
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dsvg', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials_shaded', '.svg')])



% Activity Figure with shaded error bars and limit
figure
hold on

xlabel('Time (min)')
ylabel(strcat('Normalized Average Activity (a.u.)'))
title (strcat('Normalized Average Fish Activity (All Trials) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

% x_axis = ((-length_of_trials + 2):(length_of_trials-2)); % exclude 2min from either side
x_axis = ((-length_of_trials):(2*length_of_trials)); % *** 2018a change ***


for i=1:number_of_genos
    
    if i == 1
        c = [0 0 1]; %blue
        s = 'b';
    elseif (i == 2 && number_of_genos == 2)
        c =  [1 0 0]; %red
        s = 'r';
    elseif (i == 2 && number_of_genos > 2)
        c = [0 1 1]; %cyan
        s = 'c';
    elseif (i == 3)
        c =  [1 0 0]; %red
        s = 'r';
    elseif i == 4
        c = [1 0 1]; %magenta/purple
        s = 'm';
    end
        
    sem = (nanstd(genox.all_trials_activity_normalized{i}, 0, 2)/sqrt(size(genox.all_trials_activity_normalized{i},2)))';
    shadedErrorBar(x_axis, genox.all_trials_mean_activity_normalized{i}, sem, 'lineprops', s);
    ylim([0 4])
    % legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% % axis and legend
% legend(legendname, 'Location','NorthWest')
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Average_Activity_All_Trials_shaded_Y-limit', '.fig')))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials_shaded_Y-limit', '.tif')])
% save svg
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dsvg', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Activity_All_Trials_shaded_Y-limit', '.svg')])

%% **************************************************    Sleep and Normalized Sleep    **************************************************

% Generate individual trials sleep graph data
% The difference between the graphs and the data above, is that the graphs
% includes the lights-on spike, which is not used for the analysis

genox.SideBySideSleepNormal_All_T =[];

for k=1:number_of_trials
    
    % First get the BL and PS SLEEP regions
    for i = 1:number_of_genos %loop through all genos
        genox.BLsleepT{k}{i} = genox.sleep{i}(BLstartT{k}:BLendT{k},:);
        genox.PSsleepT{k}{i} = genox.sleep{i}(PSstartT{k}:PSendT{k},:);
    end
    
    % Mean of Mean BL sleep for each geno
    for i = 1:number_of_genos %loop through all genos
        genox.BLsleepMeanOfMeanT{k}{i} = mean((sum(genox.BLsleepT{k}{i})/(BLendT{k} - BLstartT{k} +1)), 2); % Sum BL sleep for each fish, then mean that (mean of means)
    end
    
    % Mean PS sleep for each fish
    for i = 1:number_of_genos %loop through all genos
        genox.PSsleepMeanT{k}{i} = sum(genox.PSsleepT{k}{i})/(PSendT{k} - PSstartT{k} +1); % Sum AND MEAN the sleep for each fish
    end
    
    % Divide the Mean PS sleep of each fish with the mean BL sleep for the genotype
    % to get the normalized sleep
    for i = 1:number_of_genos %loop through all genos
        genox.PS_BLsleepNormalT{k}{i} = (genox.PSsleepMeanT{k}{i}/genox.BLsleepMeanOfMeanT{k}{i})';
    end
    
    genox.SideBySideSleepNormalT{k} = padcat(genox.PS_BLsleepNormalT{k}{1}, genox.PS_BLsleepNormalT{k}{2});       
    genox.SideBySideSleepNormal_All_T = cat(1, genox.SideBySideSleepNormal_All_T, genox.SideBySideSleepNormalT{k});


end

%% Sleep graph

% Get all the trial sleep slices next to each other
genox.all_trials_sleep{i} = [];
for i = 1:number_of_genos
    for k = 1:number_of_trials
        previous_slices = genox.all_trials_sleep{i};
        current_slice = genox.sleep{i}(graph_time_start{k}:graph_time_end{k}, :);
        genox.all_trials_sleep{i} = cat(2, previous_slices, current_slice);
    end
end


% Normalize with graph_BL_sleep
genox.all_trials_sleep_normalized{i} = [];
genox.all_trials_mean_sleep_normalized{i} = [];
for i = 1:number_of_genos
    genox.graph_BL_sleep{i} = mean2(genox.all_trials_sleep{i}(3:(length_of_trials - 2), :)); % 3:(length_of_trials - 2) is BLstart to BLend
    genox.all_trials_sleep_normalized{i} = genox.all_trials_sleep{i}/genox.graph_BL_sleep{i};
    genox.all_trials_mean_sleep_normalized{i} = mean(genox.all_trials_sleep{i}/genox.graph_BL_sleep{i}, 2);
end

% Make figures

% Sleep figure with no error bars
figure
hold on

xlabel('Time (min)')
ylabel(strcat('Normalized Average Sleep (a.u.)'))
title (strcat('Normalized Average Fish Sleep (All Trials) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

% *** 2018a change ***
% x_axis = ((-length_of_trials + 2):(length_of_trials-2)); % exclude 2min from either side
x_axis = ((-length_of_trials):(2 * length_of_trials)); % 

for i=1:number_of_genos
    plot(x_axis, genox.all_trials_mean_sleep_normalized{i})
    legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% axis and legend
legend(legendname, 'Location','NorthWest')
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Average_Sleep_All_Trials', '.fig')))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Sleep_All_Trials', '.tif')])
% save svg
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dsvg', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Sleep_All_Trials', '.svg')])


% Sleep Figure with shaded error bars
figure
hold on

xlabel('Time (min)')
ylabel(strcat('Normalized Average Sleep (a.u.)'))
title (strcat('Normalized Average Fish Sleep (All Trials) ', filename(1:end-9)), 'Interpreter', 'none')
set(gca,'TickDir','out',...
    'FontSize',12, 'LineWidth',2, 'FontName', 'Arial')

% *** 2018a change ***
% x_axis = ((-length_of_trials + 2):(length_of_trials-2)); % exclude 2min from either side
x_axis = ((-length_of_trials):(2 * length_of_trials)); % 


for i=1:number_of_genos
    
    if i == 1
        c = [0 0 1]; %blue
        s = 'b';
    elseif (i == 2 && number_of_genos == 2)
        c =  [1 0 0]; %red
        s = 'r';
    elseif (i == 2 && number_of_genos > 2)
        c = [0 1 1]; %cyan
        s = 'c';
    elseif (i == 3)
        c =  [1 0 0]; %red
        s = 'r';
    elseif i == 4
        c = [1 0 1]; %magenta/purple
        s = 'm';
    end
        
    sem = (nanstd(genox.all_trials_sleep_normalized{i}, 0, 2)/sqrt(size(genox.all_trials_sleep_normalized{i},2)))';
    shadedErrorBar(x_axis, genox.all_trials_mean_sleep_normalized{i}, sem, 'lineprops', s);
    % legendname(i) = cellstr(genox.name{i}); % and to update the legend
end

% % axis and legend
% legend(legendname, 'Location','NorthWest')
% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Average_Sleep_All_Trials_shaded', '.fig')))
% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Sleep_All_Trials_shaded', '.tif')])
% save svg
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 32 15])
print('-dsvg', '-r300', [working_path, '_analysis_output/', filename(1:end-9), strcat('_Opto_Normalized_Sleep_All_Trials_shaded', '.svg')])


%% Clear outliers from last 

genox.side_by_side_normalized_activity_4xIQR = outlier_elim_xIQR(genox.SideBySideActivityNormal_All_T, 3);
genox.side_by_side_normalized_sleep_4xIQR = outlier_elim_xIQR(genox.SideBySideSleepNormal_All_T,3);



%% Save 
save(strcat(working_path, '_analysis_output/',filename(1:end-9),'.mat'),'genox')

toc