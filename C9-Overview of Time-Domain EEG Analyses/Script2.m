%{
    Problem Intro.
    Loop through each electrode and find the peak time of the ERP between 100 and 400 ms. 
    Store these peak times in a separate variable and then make a topographical plot of the peak times 
    (that is, the topographical map will illustrate times in milliseconds, not activity at peak times). 
    Include a color bar in the figure and make sure to show times in milliseconds from time 0 
    (not, for example, time in milliseconds since 100 ms or indices instead of milliseconds). 
    What areas of the scalp show the earliest and the latest peak responses to the stimulus within this window?
%}

load ../data/sampleEEGdata.mat;

start_time = 100;
end_time = 400;
time_array = EEG.times;
ERP_every_nodes = mean(EEG.data, 3);
peak_time = zeros(EEG.nbchan, 1);
time_index = find(time_array >= start_time & time_array <= end_time);


% for ch_num = 1:EEG.nbchan
%     [peak_value, time] = max(ERP_every_nodes(ch_num, time_index(1):time_index(end)));
%     peak_time(ch_num) = time_array(time + time_index(1) - 1);
% end

for ch_num = 1:EEG.nbchan
    peak_value = -100000000;
    for time = time_index(1):time_index(end)
        if ERP_every_nodes(ch_num, time) > peak_value
            peak_value = ERP_every_nodes(ch_num, time);
            peak_time(ch_num) = time_array(time);
        end;
    end;
end;


topoplot(peak_time, EEG.chanlocs, 'maplimits', [min(peak_time), max(peak_time)]);
cbar('vert',0,[min(peak_time) max(peak_time)]);