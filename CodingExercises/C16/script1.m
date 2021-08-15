%{

Pick one electrode and compute a time-frequency map of power using both the multitaper method and the short-time FFT. 
Store all of the power values for all of the trials.

%}

load ../../data/sampleEEGdata.mat;

electrode2plot = 'p7';
electrode2plotidx = strcmpi(electrode2plot,{EEG.chanlocs.labels});

timewin = 400;
times2save = -500:10:1000;
times2saveidx = zeros(length(times2save),1);
for ti = 1:length(times2save)
    [junk, times2saveidx(ti)] = min(abs(EEG.times - times2save(ti)));
end

pts_in_timewin = round((timewin/1000) * EEG.srate);
hann_win = 0.5 * (1 - cos(2 * pi * (0:(pts_in_timewin - 1))/(pts_in_timewin - 1)));
freqs2save = linspace(0, EEG.srate/2, floor(pts_in_timewin/2) + 1);

stFFT_res = zeros(length(freqs2save), length(times2save), size(EEG.data, 3));

for ti = 1:length(times2saveidx)
    start_idx = times2saveidx(ti) - floor(pts_in_timewin/2);
    end_idx = times2saveidx(ti) + floor(pts_in_timewin/2) - mod(pts_in_timewin + 1, 2);

    temp_dat = squeeze(EEG.data(electrode2plotidx, start_idx:end_idx, :));
    taper_dat = bsxfun(@times, temp_dat, hann_win');
    fft_dat = fft(taper_dat, [], 1);

    stFFT_res(:, ti, :) = fft_dat(1:floor(pts_in_timewin/2) + 1, :);
end

stFFT_power = squeeze(mean(abs(stFFT_res).^2, 3));

figure
% subplot(221)
imagesc(times2save, freqs2save, zscore(stFFT_power))
colormap('jet')
set(gca, 'YDir', 'normal');
colorbar();
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Power of FFT')