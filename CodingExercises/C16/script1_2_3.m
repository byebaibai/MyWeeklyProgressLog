%{

Pick one electrode and compute a time-frequency map of power using both the multitaper method and the short-time FFT. 
Store all of the power values for all of the trials.

Next, compute a time-frequency map of signal-to-noise ratio. 
The signalto-noise ratio of power is discussed more in chapter 18, 
but it can be estimated as the average power at each time-frequency point across trials, 
divided by the standard deviation of power at each timefrequency point across trials.

Make time-frequency plots of power and signal-to-noise ratio from the two methods. 
Make another plot in which you subtract the signal-to-noise plots between the two methods. 
Are there any noticeable differences between the signal-to-noise results when the multitaper method versus the short-time FFT is used?

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

electrode2plot = 'P7';
electrode2plotidx = strcmpi(electrode2plot,{EEG.chanlocs.labels});

timewin = 400;
times2save = -300:35:1000;
times2saveidx = zeros(length(times2save));
for ti = 1:length(times2save)
    [junk, times2saveidx(ti)] = min(abs(EEG.times - times2save(ti)));
end

baseline_start = -300;
baseline_end = 0;
[junk, baseline_start_idx] = min(abs(times2save - baseline_start));
[junk, baseline_end_idx] = min(abs(times2save - baseline_end));

pts_in_timewin = round((timewin/1000) * EEG.srate);
hann_win = 0.5 * (1 - cos(2 * pi * (0:pts_in_timewin - 1)/(pts_in_timewin - 1)));

num_tapers = 5;
tapers = dpss(pts_in_timewin, num_tapers); % this line will crash without matlab signal processing toolbox

freqs2save = linspace(0, EEG.srate/2, floor(pts_in_timewin/2) + 1);

stFFT_res = zeros(length(freqs2save), length(times2save), size(EEG.data, 3));

multitaper_power = zeros(length(freqs2save), length(times2save));
multitaper_power_std = zeros(length(freqs2save), length(times2save));

for ti = 1:length(times2saveidx)
    start_idx = times2saveidx(ti) - floor(pts_in_timewin/2);
    end_idx = times2saveidx(ti) + floor(pts_in_timewin/2) - mod(pts_in_timewin + 1, 2);

    temp_dat = squeeze(EEG.data(electrode2plotidx, start_idx:end_idx, :));
    taper_dat = bsxfun(@times, temp_dat, hann_win');
    fft_dat = fft(taper_dat, [], 1);

    taperpow = zeros(floor(pts_in_timewin/2) + 1, num_tapers, size(EEG.data, 3));
    for tapi = 1:num_tapers
        multi_taper_dat = bsxfun(@times, temp_dat, tapers(:, tapi));
        multi_taper_fft = fft(multi_taper_dat, pts_in_timewin)/pts_in_timewin; % NOT UNDERSTANDED...
        multi_taper_fft = multi_taper_fft(1:floor(pts_in_timewin/2)+1, :);
        taperpow(:, tapi, :) = squeeze(taperpow(:, tapi, :)) + multi_taper_fft .* conj(multi_taper_fft);
    end
    
    stFFT_res(:, ti, :) = fft_dat(1:floor(pts_in_timewin/2) + 1, :);
    multitaper_power(:, ti) = mean(mean(taperpow, 3), 2);
    multitaper_power_std(:, ti) = std(mean(taperpow, 2), 0, 3);
end

stFFT_power = squeeze(mean(abs(stFFT_res).^2, 3));
stFFT_power_std = squeeze(std(abs(stFFT_res).^2, 0, 3));
stFFT_baseline_power = mean(stFFT_power(:, baseline_start_idx:baseline_end_idx), 2);

multitaper_baseline_power = mean(multitaper_power(:, baseline_start_idx:baseline_end_idx), 2);

stFFT_dB = 10 * log10(bsxfun(@rdivide, stFFT_power, stFFT_baseline_power));
stFFT_SNR = stFFT_power ./ stFFT_power_std;
multitaper_dB = 10 * log10(bsxfun(@rdivide, multitaper_power, multitaper_baseline_power));
multitaper_SNR = multitaper_power ./ multitaper_power_std;

figure

subplot(221)
contourf(times2save, freqs2save, multitaper_dB, 30, 'linecolor','none')
colormap('jet')
set(gca, 'YDir', 'normal', 'clim', [-3 3]);
colorbar();
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Power of multitaper')

subplot(222)
contourf(times2save, freqs2save, stFFT_dB, 30, 'linecolor','none')
colormap('jet')
set(gca, 'YDir', 'normal', 'clim', [-3 3]);
colorbar();
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Power of FFT')

subplot(223)
contourf(times2save, freqs2save, multitaper_SNR, 30, 'linecolor','none')
colormap('jet')
set(gca, 'YDir', 'normal', 'clim', [.4 1.6]);
colorbar();
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('SNR of multitaper')

subplot(224)
contourf(times2save, freqs2save, stFFT_SNR, 30, 'linecolor','none')
colormap('jet')
set(gca, 'YDir', 'normal', 'clim', [.4 1.6]);
colorbar();
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('SNR of FFT')