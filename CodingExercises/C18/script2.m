%{
Select five time points and create topographical maps of power with 
and without baseline normalization at each selected time-frequency point. 
You should have time in columns and with/without baseline normalization in rows. 
Use separate figures for each frequency. The color scaling should be the same for 
all plots over time within a frequency, but the color scaling should be different 
for with versus without baseline normalization and should also be different for 
each frequency.
%}

load ../../data/sampleEEGdata.mat;

times2plot = 0:100:400;
baseline_start_time = -300;
baseline_end_time = 0;

[junk, baseine_start_idx] = min(abs(EEG.times - baseline_start_time));
[junk, baseline_end_idx] = min(abs(EEG.times - baseline_end_time));

times2plotidx = zeros(size(times2plot));
for ti = 1:length(times2plot)
    [junk, times2plotidx(ti)] = min(abs(EEG.times - times2plot(ti)));
end
freq2plot = 8; % in Hz


% for freq = 4:1:90
% freq2plot = freq
%% Hilbert
filter_interval = 2;
data2filt = reshape(EEG.data, size(EEG.data, 1), []);
filtered_dat = reshape(eegfilt(data2filt, EEG.srate, freq2plot - filter_interval, freq2plot + filter_interval), size(EEG.data));
hilbert_res = zeros(size(filtered_dat));
for ci = 1:size(EEG.data, 1)
    hilbert_res(ci, :, :) = squeeze(filtered_dat(ci, :, :));
end

hilbert_power = squeeze(mean(abs(hilbert_res).^2, 3));
hilbert_baseline_power = squeeze(mean(hilbert_power(:, baseline_start_idx:baseline_end_idx), 2));
hilbert_db = 10 * log10(bsxfun(@rdivide, hilbert_power, hilbert_baseline_power));

figure
for ti = 1:length(times2plotidx)
    subplot(2,5,ti);
    topoplot(hilbert_power(:, times2plotidx(ti)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');
    titleStr = sprintf('Raw power, %dms', times2plot(ti));
    title(titleStr);
end

for ti = 1:length(times2plotidx)
    subplot(2,5,ti + 5);
    topoplot(hilbert_db(:, times2plotidx(ti)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');
    titleStr = sprintf('DB power, %dms', times2plot(ti));
    title(titleStr);
end
% end