%% Housecleaning
tic
clear
close all;

%% Import viewpoint data file
input_file = '241217_6A_6B_tph2ChR2_nf1ab_01_UTF8.txt';
data = readtable(strcat('../matlab_data_processed/', input_file), 'Delimiter', 'tab');

%% Set up clean_data table
% set up variables
vars = cell(1, 197);
vars{1} = 'starttimes';
vars{2} = 'endtimes';
vars{195} = 'values_are';
vars{196} = 'middur_plus_burdur';
vars{197} = 'CLOCK';
for i = 3:194
    vars{i} = strcat('FISH', num2str(i-2));
end
% set up variable types
types = cell(1,197);
for i = 1:197
    types{i} = 'double';
end

% find the number of observations
rows = height(data(strcmp(data.location, 'c1'), :));

% and now assemble the table
clean_data = table('Size', [rows, 197], 'VariableTypes', types ,'VariableNames',vars);

%% get the middur data
for i = 1:192
        current_fish = strcat('c', num2str(i,'%01.f'));
        temp_data = sortrows(data(strcmp(data.location, current_fish), :), {'start'});
        merge = table(temp_data.middur + temp_data.burdur);
        clean_data(:, i+2) = merge; % column 14 is middur, 16 is burdur
end

%% get starttimes and end times
clean_data(:, 1) = temp_data(:,7); % column 7 is starttimes
clean_data(:, 2) = temp_data(:,8); % column 8 is endtimes

for i = 1 : rows
    [h,m,s] = hms(temp_data.sttime(i));
    if h >= 9
        time_stamp = h + m/60 + s/3600 - 9; % if >= 9 subtract 9
    else
        time_stamp = h + m/60 + s/3600 + 15; % else add 15 since you are after midnight
    end
    clean_data.CLOCK(i) = time_stamp;

end
clean_data(end,:) = []; % last line of data is dirty (<60 sec)

%% save output file
output_file = strcat(input_file(1:end-4), '_DATA_midburdur.txt');
writetable(clean_data, strcat('../matlab_data_processed/', output_file), 'Delimiter', '\t');
        
  

toc
    




    
