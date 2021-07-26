%{
    Problem Intro.
    Compute the ERP at each electrode. 
    Select five time points at which to show topographical plots (e.g., 0 to 400 ms in 100-ms steps). 
    In one figure, make a series of topographical plots at these time points. 
    To increase the signal-to-noise ratio, 
    make each plot show the average of activity from 20 ms before until 20 ms after each time point. 
    For example, the topographical plot from 200 ms should show average activity from 180 ms until 220 ms. 
    Indicate the center time point in a title on each subplot.
%}

load ../data/sampleEEGdata.mat;

ERP_every_nodes = mean(EEG.data, 3);
time_array = EEG.times;
time_offset = 20;
max_abs_value = -1000000;

tiledlayout(1, 5);
for index = linspace(0, 400, 5)
    nexttile;
    base_index = index + 1;
    time_range = find(time_array >= index - time_offset & time_array <= index + time_offset);
    ERP_one_plot = mean(ERP_every_nodes(:, time_range(1) : time_range(end)), 2);
    topoplot(ERP_one_plot, EEG.chanlocs, 'style', 'map', 'electrodes', 'off', 'maplimits', [-10, 10]);
    max_abs_value = max(max_abs_value, abs(ERP_one_plot));

    title_s = sprintf("ERP at %d ms", index);
    set(get(gca, 'Title'), 'String', title_s);    
end

cbar('vert',0,[-10 10]);