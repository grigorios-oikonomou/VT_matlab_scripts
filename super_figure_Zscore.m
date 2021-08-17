function super_figure_Zscore(day_number, genox)

%% Some default settings
set(0, 'DefaultFigureRenderer', 'painters');
Zsummary_figure = figure('pos',[500 500 1600 1200]);


%% Activity 60 min
ax1 = subplot(4,6,[1 2]);
xlabel('Zeitgeber Time (h)')
ylabel('Average Activity (s/h)')
title(strcat('Activity 1h bins -'))
title(strcat('Activity 1h bins -', genox.experiment)) % to reset
fish_plot_shaded(genox.onehourtime, genox.onehour, genox) % fish_plot(time, data, rest)


%% Activity 1 min
ax19 = subplot(4,6,[19 20]);
hold on
xlabel('Zeitgeber Time (hours)')
ylabel('Average Activity (s/min)')
title(strcat('Activity 1min bins -', genox.experiment))
fish_plot(genox.OneMinuteTime, genox.OneMinute, genox) % fish_plot(time, data, rest)

%% Waking activity 60 min
ax2 = subplot(4,6,[3 4]);
hold on
xlabel('Zeitgeber Time (hours)')
ylabel('Waking Activity (s/awake min)')
title(strcat('Waking Activity | 1 h bins ',genox.data_filename(1:end-9))) 
fish_plot_shaded(genox.onehourtime, genox.avewakechart1h, genox) % fish_plot(time, data, rest)

%% Waking Activity 1 min
ax21 = subplot(4,6,[21 22]);
hold on
xlabel('Zeitgeber Time (hours)')
ylabel('Average Waking Activity (s/awake min')
title(strcat('Waking Activity | 1 min bins ',genox.data_filename(1:end-9)))
fish_plot(genox.OneMinuteTime, genox.avewaking, genox) % fish_plot(time, data, rest)

%% Sleep graph 60 min
ax3 = subplot(4,6,[5 6]);
hold on
xlabel('Zeitgeber Time (h)')
ylabel('Sleep (min/h)')
title (strcat('Sleep | 60 min bins ', genox.data_filename(1:end-9)))
fish_plot_shaded(genox.onehourtime, genox.sleepchart1h, genox) % fish_plot(time, data, rest)
legend('Location','SouthEast');

%% Sleep 1 min
ax23 = subplot(4,6,[23 24]);
hold on
xlabel('Zeitgeber Time (h)')
ylabel('Sleep (min/min)')
title (strcat('Sleep | 1 min bins ', genox.data_filename(1:end-9)))
fish_plot(genox.OneMinuteTime, genox.sleep, genox) % fish_plot(time, data, rest)
legend('Location','SouthEast');

%% Boxplots with extreme outliers not shown

ax7 = subplot(4,6,7);
boxplot(genox.Zsummary.hourly_activity.day{day_number},'Notch','on','width', 0.75) 
title(strcat("Z-score Hourly Activity Day ", num2str(day_number), ' (s/h)'))
var = genox.Zsummary.hourly_activity.day{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax7.YLim = [minY, maxY];

ax8 = subplot(4,6,8);
boxplot(genox.Zsummary.waking_activity.day{day_number},'Notch','on','width', 0.75) 
title(strcat("Z-score Waking Activity Day ", num2str(day_number), ' s/awake min'))
var = genox.Zsummary.waking_activity.day{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax8.YLim = [minY, maxY];

ax9 = subplot(4,6,9);
boxplot(genox.Zsummary.hourly_sleep.day{day_number},'Notch','on','width', 0.75) 
title(strcat("Z-score Hourly Sleep Day ", num2str(day_number), ' (min/h)'))
var = genox.Zsummary.hourly_sleep.day{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax9.YLim = [minY, maxY];

ax10 = subplot(4,6,10);
boxplot(genox.Zsummary.hourly_sleep_bout_number.day{day_number},'Notch','on','width', 0.75) 
title(strcat("Z-score Sleep Bouts Day ", num2str(day_number), ' (per h)'))
var = genox.Zsummary.hourly_sleep_bout_number.day{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax10.YLim = [minY, maxY];

ax11 = subplot(4,6,11);
boxplot(genox.Zsummary.average_sleep_bout_length.day{day_number},'Notch','on','width', 0.75) 
title(strcat("Z-score Sleep Bouts Length Day ", num2str(day_number), ' (min)'))
var = genox.Zsummary.average_sleep_bout_length.day{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax11.YLim = [minY, maxY];


ax13 = subplot(4,6,13);
boxplot(genox.Zsummary.hourly_activity.night{day_number},'Notch','on','width', 0.75) 
ax13.Color = [0.9 0.9 0.9];
title(strcat("Z-score Hourly Activity Night ", num2str(day_number), ' (s/h)'))
var = genox.Zsummary.hourly_activity.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax13.YLim = [minY, maxY];

ax14 = subplot(4,6,14);
boxplot(genox.Zsummary.waking_activity.night{day_number},'Notch','on','width', 0.75) 
ax14.Color = [0.9 0.9 0.9];
title(strcat("Z-score Waking Activity Day Night ", num2str(day_number), ' s/awake min'))
var = genox.Zsummary.waking_activity.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax14.YLim = [minY, maxY];

ax15 = subplot(4,6,15);
boxplot(genox.Zsummary.hourly_sleep.night{day_number},'Notch','on','width', 0.75) 
ax15.Color = [0.9 0.9 0.9];
title(strcat("Z-score Hourly Sleep Night ", num2str(day_number), ' (min/h)'))
var = genox.Zsummary.hourly_sleep.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax15.YLim = [minY, maxY];

ax16 = subplot(4,6,16);
boxplot(genox.Zsummary.hourly_sleep_bout_number.night{day_number},'Notch','on','width', 0.75) 
ax16.Color = [0.9 0.9 0.9];
title(strcat("Z-score Sleep Bouts Night ", num2str(day_number), ' (per h)'))
var = genox.Zsummary.hourly_sleep_bout_number.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax16.YLim = [minY, maxY];

ax17 = subplot(4,6,17);
boxplot(genox.Zsummary.average_sleep_bout_length.night{day_number},'Notch','on','width', 0.75) 
ax17.Color = [0.9 0.9 0.9];
title(strcat("Z-score Sleep Bouts Length Night ", num2str(day_number), ' (min)'))
var = genox.Zsummary.average_sleep_bout_length.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax17.YLim = [minY, maxY];

ax18 = subplot(4,6,18);
boxplot(genox.Zsummary.sleep_latency.night{day_number},'Notch','on','width', 0.75) 
ax18.Color = [0.9 0.9 0.9];
title(strcat("Z-score Latency Night ", num2str(day_number), ' (min)'))
var = genox.Zsummary.sleep_latency.night{day_number};
p75 = prctile(var,75);
p25 = prctile(var,25);
IQR = p75 - p25;
tops = p75 + 2*IQR;
maxY = max(tops);
bottoms = p25 - 2*IQR;
minY = min(bottoms);
ax18.YLim = [minY, maxY];

end








