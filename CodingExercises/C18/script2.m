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

[junk, baseline_start_idx] = min(abs(EEG.times - baseline_start_time));
[junk, baseline_end_idx] = min(abs(EEG.times - baseline_end_time));

times2plotidx = zeros(size(times2plot));
for ti = 1:length(times2plot)
    [junk, times2plotidx(ti)] = min(abs(EEG.times - times2plot(ti)));
end
freq2plot = 6; % in Hz

% Hilbert Filtering
% filter_interval = 2;
% data2filt = reshape(EEG.data, size(EEG.data, 1), []);
% filtered_dat = reshape(eegfilt(data2filt, EEG.srate, freq2plot - filter_interval, freq2plot + filter_interval), size(EEG.data));
% hilbert_res = zeros(size(filtered_dat));
% for ci = 1:size(EEG.data, 1)
%     hilbert_res(ci, :, :) = squeeze(hilbert(filtered_dat(ci, :, :)));
% end

% hilbert_power = squeeze(mean(abs(hilbert_res).^2, 3));
% hilbert_baseline_power = squeeze(mean(hilbert_power(:, baseline_start_idx:baseline_end_idx), 2));
% hilbert_db = 10 * log10(bsxfun(@rdivide, hilbert_power, hilbert_baseline_power));

%% Wavelet
t = 1.23;
times = linspace(-t, t, size(EEG.data, 2));
sin_wave = exp(2 * pi * 1i * freq2plot .* times);
std_gauss = 5/(2*pi*freq2plot);
exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
wavelet = sin_wave .* exp_wave;
wavelet_res = zeros(size(EEG.data));
for eli = 1:size(EEG.data, 1)
    for tri = 1:99
        wavelet_res(eli, :, tri) = conv(EEG.data(eli, :, tri), wavelet, 'same');
    end
end
wavelet_power = squeeze(mean(abs(wavelet_res).^2, 3));
wavelet_baseline_power = squeeze(mean(wavelet_power(:, baseline_start_idx:baseline_end_idx), 2));
wavelet_db = 10 * log10(bsxfun(@rdivide, wavelet_power, wavelet_baseline_power));

figure
for ti = 1:length(times2plotidx)
    subplot(2,5,ti);
    % not sure about this scaling, Mike said color scaling for raw power is [-120 120], but I don't know how to do this
    topoplot(wavelet_power(:, times2plotidx(ti)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off', 'maplimits',[-26000 26000]);
    titleStr = sprintf('Raw power, %dms', times2plot(ti));
    title(titleStr);
end

for ti = 1:length(times2plotidx)
    subplot(2,5,ti + 5);
    topoplot(wavelet_db(:, times2plotidx(ti)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off', 'maplimits',[-3 3]);
    titleStr = sprintf('DB power, %dms', times2plot(ti));
    title(titleStr);
end
