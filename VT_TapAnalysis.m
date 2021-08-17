%% Clean up
close all
clear all

%% Define the path you will be working in:
tic
working_path = '/Users/grigoriosoikonomou/Documents/Prober_Lab/Videotrackers/analysis_software/analysis2/';
genox.ShamWindow = 5; % this is how many sec before real tap you measure the sham tap


%% Import genotype list(listname) and dataset

listname = '170202_5A_genotype_3.txt';
% filename = '170202_5A_DATA.txt';

filename = '170202_5A_5B_big_test_DATA_middur_plus_burdur_clean_old_reformat.txt';

gridmap = importdata(strcat(working_path,'_genotype_files/', listname), '\t', 2);
dataset = importdata(strcat(working_path,'_data_files/', filename), '\t', 2);

%% Import the names of the genotypes and sort the data in different cells based on genotypes 

number_of_genos = str2double(listname(20));
% grab the 20th character from the name of the genotype file (listname) and
% then convert that to a double; this way I have the number of genotypes
% (number_of_genos)

for i=1:number_of_genos
    genox.name{i} = gridmap.colheaders{i};
    genox.data{i} = dataset.data(:, gridmap.data(~isnan(gridmap.data(:,i)),i)+2);
end

%% Save the fish numbers!!!

for i=1:number_of_genos
    genox.fishID{i} = gridmap.data(~isnan(gridmap.data(:,i)), i);
end

%% Set the timestamp

genox.timestamp = dataset.data(:,97); % well of LED +2 = number you want here
%genox.timestamp(find(genox.timestamp ~= 1)) = 0
genox.timestamp(genox.timestamp == 1) = 1;
genox.timestamp(genox.timestamp == 2) = 1;
%% Get the zeitgeber info from the dataset file (last column)

genox.zeitgeber = dataset.data(:, end);

%% Get the light schedule info from the zeitgeber info

% lights are on at time 0, lights are off at time 14

genox.lightschedule = genox.zeitgeber; % first just copy the zeitgeber

genox.lightschedule(genox.lightschedule < 14) = 0; % if <14 then lights are off, therefore 0;
genox.lightschedule(genox.lightschedule >= 14) = 1; % if >= 14 lights are on, therefore 1;

% Remember: the stamp is 0 for day, and 1 for night.


%% Get the genox.time in hours

% count the seconds:
for j=1:length(genox.lightschedule)
    genox.timesec(j) = (j);
end


genox.time = genox.timesec/3600; % now convert the seconds to hours


%% Binary movement

% convert all activity to 0 for no movement, 1 for movement
% genox.active{i} is an array of the original data, split by genotypes, with
% all activity above 0.1 converted to 1 (active) and all activitity below
% 0.1 converted to 0 (inactive)
for i = 1:number_of_genos
    genox.active{i} = genox.data{i}; % first just copy the data
    genox.active{i}(genox.data{i} <= 0.1 ) = 0; %0.1 is the minimum value given by the VT; if equal or below that, no movement;
    genox.active{i}(genox.data{i} > 0.1 ) = 1; % if above that, movement;
end


% genox.activetotal{i} is a vector (1 column) with the activity status
% (0/1) of all fish for each timepoint added 
% (essentially the genox.active array collapsed onto the first collumn)

for i = 1:number_of_genos %loop through all genos
    genox.activetotal{i} = sum(genox.active{i}')';
end


%% Graph #1: The fish activity graph

figure
hold on
xlabel('Time since start')
ylabel('Number of fish moving')
title (strcat('Fish reponse ', filename(1:end-9)), 'Interpreter', 'none')


for j = 1:number_of_genos

% First to define the colors that will be used for the different genotypes:
% copied from JR
 
    if j == 1
        c = [0 0 1];
    end
    if j == 2
        c = [1 0 0];
    end
    if j == 3 
        c = [1 1 0];
    end
    if j == 4 
        c = [1 0 1];
    end
    if j == 5
        c = [0 1 0];
    end
    if j == 6 
        c = [0 1 1];
    end
    if j == 7
        c = [0 0 0];
    end
    if j == 8 
        c = [0.5 0.5 0.5];
    end
    if j == 9
        c = [.5 .2 .8];
    end 
    if j == 10 
        c = [0.5 .2 .2];
    end
    if j > 10
        c = [1 1 1];
    end
    
    % now to plot each genotype (j)
    plot(genox.time, genox.activetotal{j}',...
        'color', c,...
        'marker','.')
    


    % and to update the legend
    legendname(j) = cellstr(genox.name{j});
end

% axis and legend
axis([genox.time(1) genox.time(end) 0 96])
legend(legendname)

% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), '_Activity_Plot.fig'))

% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), '_Activity_Plot.tif'])



%% Filtering through the timestamp


for i = 1:number_of_genos
        genox.activeclean{i} = genox.activetotal {i}; % first just coppy the data
        genox.activeclean{i} =  genox.activeclean{i}.*genox.timestamp; % then multiply with the timestamp vector (1 for active, 0 for incative) 
        
end



%% Graph #2: The clean active graph

figure
hold on
xlabel('Time since start')
ylabel('Number of fish moving above threshold')
title (strcat('Fish reponse ', filename(1:end-9)), 'Interpreter', 'none')


for j = 1:number_of_genos

% First to define the colors that will be used for the different genotypes:
% copied from JR
 
    if j == 1
        c = [0 0 1];
    end
    if j == 2
        c = [1 0 0];
    end
    if j == 3 
        c = [1 1 0];
    end
    if j == 4 
        c = [1 0 1];
    end
    if j == 5
        c = [0 1 0];
    end
    if j == 6 
        c = [0 1 1];
    end
    if j == 7
        c = [0 0 0];
    end
    if j == 8 
        c = [0.5 0.5 0.5];
    end
    if j == 9
        c = [.5 .2 .8];
    end 
    if j == 10 
        c = [0.5 .2 .2];
    end
    if j > 10
        c = [1 1 1];
    end
    
    % now to plot each genotype (j)
    plot(genox.time, genox.activeclean{j}',...
        'color', c,...
        'marker','.')
    


    % and to update the legend
    legendname(j) = cellstr(genox.name{j});
end

% axis and legend
axis([genox.time(1) genox.time(end) 0 96])
legend(legendname)

% save the fig
hgsave(strcat(working_path, '_analysis_output/', filename(1:end-9), '_CleanActive_Plot.fig'))

% save the tif
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 15 15])
print('-dtiff', '-r300', [working_path, '_analysis_output/', filename(1:end-9), '_CleanActive_Plot.tif'])




%% Finding the tap times from the timestamp! Also the Sham tap times

for i = 1:number_of_genos
    genox.taptimesALL{i}=find(genox.timestamp>0);
end


% Making the first of a double tap into a NaN
for i = 1:number_of_genos
    genox.taptimesNaNDoubles{i} = genox.taptimesALL{i}; % first just copy the  genox.taptimesALL{i} (doubles included) to genox.taptimesNaNDoubles{i}
    for j=1:length(genox.taptimesALL{i}) - 1 % for # of ALL tap times -1
        if (genox.taptimesALL{i}(j+1) - genox.taptimesALL{i}(j)) == 1 % if subtracting j from j+1 gives you 1 (only 1 sec of difference between two tap events)
           genox.taptimesNaNDoubles{i}(j,:)= [NaN]; % then turn j into a series of NaNs
        end
    end
end

% then droping all the NaNs
for i=1:number_of_genos % for all genos
A =genox.taptimesNaNDoubles{i};
genox.taptimesNONDoubles{i} = A(all(~isnan(A),2),:); % copy everything that is not NaN
end

% then copy the genox.taptimesNaNDoubles to genox.taptimes and work with that
for i = 1:number_of_genos
    genox.taptimes{i} = genox.taptimesNONDoubles{i};
end


%determine all the Sham tap times (5sec before real tap)
for i = 1:number_of_genos
    genox.taptimes_Sham{i}= genox.taptimes{i} - genox.ShamWindow; % sham tapping is 5sec before real tapping
end
%% Determining which fish were active in the 60 sec before the tap (awake fish): genox.pretapactive

for i = 1:number_of_genos % for all genos
    for j = 1:length(genox.data{i}(1,1:end)) % for each fish/column in that geno
        for q = 1:length(genox.taptimes{i}) % for each row, that is each tap time
            genox.pretapactivity{i}(q,j) = sum(genox.data{i}(genox.taptimes{i}(q)-61:genox.taptimes{i}(q)-1,j)); % this sums 60sec rec from before tap
            
            genox.pretapactive{i} = genox.pretapactivity{i}; % first just copy the data
            genox.pretapactive{i}(genox.pretapactivity{i} <= 0.1 ) = 0; %0.1 is the minimum value given by the VT; if equal or below that, no movement
            genox.pretapactive{i}(genox.pretapactivity{i} > 0.1 ) = 1; % if above that, movement
        
        end
    end
end

%% Determining which fish were active in the 60 sec before the Sham tap (awake fish): genox.pretapactive_Sham

for i = 1:number_of_genos % for all genos
    for j = 1:length(genox.data{i}(1,1:end)) % for each fish/column in that geno
        for q = 1:length(genox.taptimes_Sham{i}) % for each row, that is each tap time
            genox.pretapactivity_Sham{i}(q,j) = sum(genox.data{i}(genox.taptimes_Sham{i}(q)-61:genox.taptimes_Sham{i}(q)-1,j)); % this sums 60sec rec from before tap
            
            genox.pretapactive_Sham{i} = genox.pretapactivity_Sham{i}; % first just copy the data
            genox.pretapactive_Sham{i}(genox.pretapactivity_Sham{i} <= 0.1 ) = 0; %0.1 is the minimum value given by the VT; if equal or below that, no movement
            genox.pretapactive_Sham{i}(genox.pretapactivity_Sham{i} > 0.1 ) = 1; % if above that, movement
        
        end
    end
end


%% Use the genox.taptimes to selectively copy from the genox.active date set
% thus creating the genox.tapactive: a subest of the data that has the
% activity only at the tapping times as opposed to all times

for i = 1:number_of_genos % for all genos
    for j = 1:length(genox.active{i}(1,1:end)) % for each fish/column in that geno
        for q = 1:length(genox.taptimes{i}) % for each row, that is each tap time
            genox.tapactive{i}(q,j) = genox.active{i}(genox.taptimes{i}(q),j); % the genox.tapactive entry is a copy of the genox.active row, 
                                                                              % guided by the genox.taptimes
        end
    end
end

%% Use the genox.taptimes_Sham to selectively copy from the genox.active date set
% thus creating the genox.tapactive_Sham: a subest of the data that has the
% activity only at the Sham tapping times as opposed to all times

for i = 1:number_of_genos % for all genos
    for j = 1:length(genox.active{i}(1,1:end)) % for each fish/column in that geno
        for q = 1:length(genox.taptimes_Sham{i}) % for each row, that is each Sham tap time
            genox.tapactive_Sham{i}(q,j) = genox.active{i}(genox.taptimes_Sham{i}(q),j); % the genox.tapactive entry is a copy of the genox.active row, 
                                                                              % guided by the genox.taptimes
        end
    end
end



%% Copy (convert) genox.pretapactive to genox.pretapwake; also for Sham


for i=1:number_of_genos
    genox.pretapwake{i}=genox.pretapactive{i};
end


for i=1:number_of_genos
    genox.pretapwake_Sham{i}=genox.pretapactive_Sham{i};
end

%% Table for each experiment and each Sham experiment

for i=1:number_of_genos
    
    %real taps
    
    genox.wake_active{i} = genox.pretapwake{i}./genox.tapactive{i};
    % wake/move
    % genox.TwakeTmove = wake and moving after tap: 1/1 = 1
    % genox.TwakeFmove = wake and not moving after tap: 1/0 = Inf
    % genox.FwakeTmove = asleep and moving after tap: 0/1 = 0
    % genox.FwakeFmove = asleep and not moving after tap: 0/0 = nan
    
    % get genox.TwakeTmove = wake and moving after tap: 1/1 = 1
    genox.TwakeTmove{i} = genox.wake_active{i}./genox.wake_active{i};
    genox.TwakeTmove{i}(isnan(genox.TwakeTmove{i}))=0;
    genox.TwakeTmoveTotal{i} = sum(genox.TwakeTmove{i}')';
    
    % get genox.TwakeFmove = wake and not moving after tap: 1/0 = Inf
    genox.TwakeFmove{i}= genox.wake_active{i};
    genox.TwakeFmove{i}(isnan(genox.TwakeFmove{i}))=0; % convert all "nan" to "0"
    genox.TwakeFmove{i}(isfinite(genox.TwakeFmove{i}))=0; % convert all numbers to "0"
    genox.TwakeFmove{i}(isinf(genox.TwakeFmove{i}))=1; % convert all "inf" to "1"
    genox.TwakeFmoveTotal{i} = sum(genox.TwakeFmove{i}')';
    
    % get genox.FwakeTmove = asleep and moving after tap: 0/1 = 0
    genox.FwakeTmove{i} = genox.wake_active{i} ;
    genox.FwakeTmove{i}(genox.FwakeTmove{i} == 1)= NaN; % convert all the "1" to "nan"
    genox.FwakeTmove{i}(isinf(genox.FwakeTmove{i}))= NaN; % convert all "inf" to "nan"
    genox.FwakeTmove{i}(isfinite(genox.FwakeTmove{i}))= 1; % convert all numbers (now only 0 left) to "1"
    genox.FwakeTmove{i}(isnan(genox.FwakeTmove{i}))= 0; % convert all "nan" to "0"
    genox.FwakeTmoveTotal{i} = sum(genox.FwakeTmove{i}')';
    
    % genox.FwakeFmove = asleep and not moving after tap: 0/0 = nan
    genox.FwakeFmove{i} = genox.wake_active{i} ;
    genox.FwakeFmove{i}(isinf(genox.FwakeFmove{i}))= 0; % convert all "inf" to "0"
    genox.FwakeFmove{i}(genox.FwakeFmove{i} == 1)= 0; % convert all "1" to "0"
    genox.FwakeFmove{i}(isnan(genox.FwakeFmove{i}))=1; % convert all "nan" to "1"
    genox.FwakeFmoveTotal{i} = sum(genox.FwakeFmove{i}')';
    
    % sham taps
    
    genox.wake_active_Sham{i} = genox.pretapwake_Sham{i}./genox.tapactive_Sham{i};
    % wake/move
    % genox.TwakeTmove = wake and moving after tap: 1/1 = 1
    % genox.TwakeFmove = wake and not moving after tap: 1/0 = Inf
    % genox.FwakeTmove = asleep and moving after tap: 0/1 = 0
    % genox.FwakeFmove = asleep and not moving after tap: 0/0 = nan
    
    % get genox.TwakeTmove = wake and moving after tap: 1/1 = 1
    genox.TwakeTmove_Sham{i} = genox.wake_active_Sham{i}./genox.wake_active_Sham{i};
    genox.TwakeTmove_Sham{i}(isnan(genox.TwakeTmove_Sham{i}))=0;
    genox.TwakeTmoveTotal_Sham{i} = sum(genox.TwakeTmove_Sham{i}')';
    
    % get genox.TwakeFmove = wake and not moving after tap: 1/0 = Inf
    genox.TwakeFmove_Sham{i}= genox.wake_active_Sham{i};
    genox.TwakeFmove_Sham{i}(isnan(genox.TwakeFmove_Sham{i}))=0; % convert all "nan" to "0"
    genox.TwakeFmove_Sham{i}(isfinite(genox.TwakeFmove_Sham{i}))=0; % convert all numbers to "0"
    genox.TwakeFmove_Sham{i}(isinf(genox.TwakeFmove_Sham{i}))=1; % convert all "inf" to "1"
    genox.TwakeFmoveTotal_Sham{i} = sum(genox.TwakeFmove_Sham{i}')';
    
    % get genox.FwakeTmove = asleep and moving after tap: 0/1 = 0
    genox.FwakeTmove_Sham{i} = genox.wake_active_Sham{i} ;
    genox.FwakeTmove_Sham{i}(genox.FwakeTmove_Sham{i} == 1)= NaN; % convert all the "1" to "nan"
    genox.FwakeTmove_Sham{i}(isinf(genox.FwakeTmove_Sham{i}))= NaN; % convert all "inf" to "nan"
    genox.FwakeTmove_Sham{i}(isfinite(genox.FwakeTmove_Sham{i}))= 1; % convert all numbers (now only 0 left) to "1"
    genox.FwakeTmove_Sham{i}(isnan(genox.FwakeTmove_Sham{i}))= 0; % convert all "nan" to "0"
    genox.FwakeTmoveTotal_Sham{i} = sum(genox.FwakeTmove_Sham{i}')';
    
    % genox.FwakeFmove = asleep and not moving after tap: 0/0 = nan
    genox.FwakeFmove_Sham{i} = genox.wake_active_Sham{i} ;
    genox.FwakeFmove_Sham{i}(isinf(genox.FwakeFmove_Sham{i}))= 0; % convert all "inf" to "0"
    genox.FwakeFmove_Sham{i}(genox.FwakeFmove_Sham{i} == 1)= 0; % convert all "1" to "0"
    genox.FwakeFmove_Sham{i}(isnan(genox.FwakeFmove_Sham{i}))= 1; % convert all "nan" to "1"
    genox.FwakeFmoveTotal_Sham{i} = sum(genox.FwakeFmove_Sham{i}')';
    
    
    %summary tables
    
    genox.summary{i} = [genox.TwakeTmoveTotal{i}, genox.TwakeFmoveTotal{i}, genox.FwakeTmoveTotal{i}, genox.FwakeFmoveTotal{i}];
    genox.summary_Sham{i} = [genox.TwakeTmoveTotal_Sham{i}, genox.TwakeFmoveTotal_Sham{i}, genox.FwakeTmoveTotal_Sham{i}, genox.FwakeFmoveTotal_Sham{i}];
    
    genox.summaryTime{i} = [genox.taptimes{i}, genox.TwakeTmoveTotal{i}, genox.TwakeFmoveTotal{i}, genox.FwakeTmoveTotal{i}, genox.FwakeFmoveTotal{i}];
    genox.summaryTime_Sham{i} = [genox.taptimes_Sham{i}, genox.TwakeTmoveTotal_Sham{i}, genox.TwakeFmoveTotal_Sham{i}, genox.FwakeTmoveTotal_Sham{i}, genox.FwakeFmoveTotal_Sham{i}];
    
    genox.ShamMovedRatio{i} = (genox.TwakeTmoveTotal_Sham{i} + genox.FwakeTmoveTotal_Sham{i})./length(genox.active{i}(1,1:end));
    genox.summaryTime_ShamMovedRatio{i} = [genox.taptimes{i}, genox.TwakeTmoveTotal{i}, genox.TwakeFmoveTotal{i}, genox.FwakeTmoveTotal{i}, genox.FwakeFmoveTotal{i}, genox.ShamMovedRatio{i}];
    genox.summaryTime_Exp_Sham{i} = [genox.summaryTime{i}, genox.summary_Sham{i}];
end



save(strcat(working_path, '_analysis_output/',filename(1:end-9),'.mat'),'genox')


toc