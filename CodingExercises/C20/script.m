%{

Pick two frequencies and compute total and non-phase-locked power from each 
electrode over time at these two frequencies. Pick two time windows, one early and one late, 
of several hundreds of milliseconds each (e.g., 100-300 ms and 500-800 ms) 
and show topographical maps of total power, nonphase-locked power, and phase-locked power 
from the average of all time points within these windows. Are there striking topographical 
differences among these results? If so, are the differences bigger or smaller in the early or the 
late time window? Why might this be the case?

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

freqs2plot = [5 20]; % Hz
timewin2plot = [100 300; 600 900];
timewin2plotidx = zeros(size(timewin2plot));
for ri = 1:size(timewin2plot, 1)
    for ci = 1:size(timewin2plot, 2)
        [~, timewin2plotidx(ri, ci)] = min(abs(EEG.times - timewin2plot(ri, ci)));
    end
end
baseline_time = [-400 -100];
[~, baseidx(1)] = min(abs(EEG.times - baseline_time(1)));
[~, baseidx(2)] = min(abs(EEG.times - baseline_time(2)));


% FFT parameters
times = -1:1/EEG.srate:1;
n_wavelet = length(times);
n_data = EEG.pnts * EEG.trials;
n_conv(1:2) = n_wavelet + n_data - 1;
n_conv(3) = n_wavelet + EEG.pnts - 1;
half_wavelet_size = (length(times) - 1)/2;

% compute ERP
erp = mean(EEG.data, 3);

% compute induced power(non-phased-locked) by substracting ERP from each trial
induced_EEG = EEG.data - erp;


% FFT of data
fft_total = fft(reshape(EEG.data, EEG.nbchan, EEG.pnts * EEG.trials), n_conv(1), 2); % total
fft_induced = fft(reshape(induced_EEG, EEG.nbchan, EEG.pnts * EEG.trials), n_conv(2),2); % non-phased-locked
fft_evoked = fft(erp, n_conv(3), 2); % erp

fft_EEG{1} = fft_total;
fft_EEG{2} = fft_induced;
fft_EEG{3} = fft_evoked;

tf_res = zeros(4, length(freqs2plot), EEG.nbchan, EEG.pnts);

for fi = 1:length(freqs2plot)
    freq2plot = freqs2plot(fi);

    wavelet = exp(2 * 1i * pi * freq2plot .* times) .* exp(-times.^2 ./ (2 * (4/(2 * pi * freq2plot))^2)) / freq2plot;
    for si=1:3
        fft_wavelet = fft(wavelet, n_conv(si));
        for chi=1:EEG.nbchan
            fft_eeg = fft_EEG{si}(chi, :);
            conv_result = ifft(fft_eeg .* fft_wavelet, n_conv(si));
            conv_result = conv_result(half_wavelet_size+1:end-half_wavelet_size);

            if si < 3
                conv_result = reshape(conv_result, EEG.pnts, EEG.trials);
                tf_res(si, fi, chi, :) = mean(abs(conv_result).^2, 2);
            else
                tf_res(si, fi, chi, :) = abs(conv_result).^2;
            end

            % db correct power
            tf_res(si, fi, chi, :) = 10*log10( squeeze(tf_res(si, fi, chi, :)) ./ mean(tf_res(si, fi, chi, baseidx(1):baseidx(2)),4) );
        end    
    end
end
tf_res(4, :, :, :) = tf_res(1, :, :, :) - tf_res(2, :, :, :);

analysis_labels = {'Total'; 'Non-phase-locked'; 'Phase-locked'};
label_index = [1, 2, 4];

figure
for fi = 1:length(freqs2plot)
    for ti = 1:length(timewin2plot)
        for li = 1:3
            data2plot = squeeze(mean(tf_res(label_index(li), fi, :, timewin2plotidx(ti, 1):timewin2plotidx(ti, 2)), 4));
            figure_idx = li + (fi - 1) * length(timewin2plot) * length(analysis_labels) + (ti - 1) * length(analysis_labels);
            subplot(length(timewin2plot)*length(freqs2plot), length(analysis_labels), figure_idx)
            topoplot(data2plot, EEG.chanlocs, 'style', 'map', 'electrodes', 'off', 'maplimits', [-2 2]);
            titleStr = sprintf('%s (%dHz, %d-%dms)',analysis_labels{li}, freqs2plot(fi), timewin2plot(ti, 1), timewin2plot(ti, 2));
            title(titleStr)
        end
    end
end

