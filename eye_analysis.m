%% 

clear all;clc;

%%Important directories

expDir = [pwd filesep];
outDir = [pwd filesep];

%% Find files

% select all subject files in the current folder
subjects = dir('subj*.txt');
subjects = {subjects(:).name};

%% loop thru the rows of raw data, col1=timestamp, col2=x, col3=y
for i = 1:length(subjects)
    raw_data = importdata(subjects{i});
    distance = zeros(length(raw_data),2);
    for j = 1:length(raw_data)
        if j ~= length(raw_data) && raw_data(j,2) ~=0 % if j is neither the last element nor 0
            x2 = raw_data(j+1,2);
            x1 = raw_data(j,2);
            y2 = raw_data(j+1,3);
            y1 = raw_data(j,3);
            distance(j,2) = sqrt((x2-x1)^2 + (y2-y1)^2); % calculate the straight line (distance) between two points on a 2-dimensional plane
        else
            distance(j,2) = 0;
        end
    end
    distance(:,1) = raw_data(:,1); %copy the timestamp from raw_data to distance
    
    first = raw_data(1,1);
    last = raw_data(length(raw_data),1);
    ref_timestamp = [first:1:last]';
    short_timestamp = distance(:,1);
    
    idx = setdiff(ref_timestamp,short_timestamp); %find the difference elements between two arrays
    missing_timestamp = zeros(length(idx),2);
    missing_timestamp(:,1) = idx; %this array contains all the missing timestamps
    combined_timestamp = [missing_timestamp;distance]; %combine two arrays
    final_timestamp = sortrows(combined_timestamp,1); %sort according to timestamps
    n = 2000; %data size in one TR (bin)
    s1 = length(final_timestamp(:,2)); %total data points
    bins = [round(s1/n)-1]; 
    y = s1 - bins.*n;
    even_bins = bins.*n;
    reshapeBins = reshape(final_timestamp(1:even_bins,2),n,bins); %reshape to N bins, each with 2000 rows
    last_bin = sum(final_timestamp(y:end,2)); %the remaining data as the last data point
    sumbins = sum(reshapeBins); %sump up
    results(:,i) = [sumbins last_bin]; %include the last datapoint and save it to matrix results
end

