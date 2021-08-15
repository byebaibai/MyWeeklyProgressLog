%{

Select one electrode and one frequency and compute power over 
time at that electrode and that frequency using complex wavelet convolution, 
filterHilbert, and the short-time FFT. 
Plot the results of these three timefrequency decomposition methods in
 different subplots of one figure. 
 Note that the scaling might be different 
 because no baseline normalization has been applied. 
 How visually similar are the results from these three methods? 
 If the results from the three methods are different, 
 how are they different, and what parameters do you think you 
 could change in the three methods to make the results look more or less similar?

%}


% Same with Solution

load ../../data/sampleEEGdata.mat;

timewin = 400; % in ms, for short-time FFT
t_s = find(EEG.times >= -200, 1);
t_e = find(EEG.times <= 800, 1, 'last');
time_length = t_e - t_s + 1;
times2save = linspace(-200, 800, time_length); % in ms
channel2plot = 47;
freq2plot = 6; % in Hz
channel2useidx = (linspace(1, 64, 64) == channel2plot);

%% Short-time FFT
times2saveidx = zeros(size(times2save));
for i = 1:length(times2save)
    [junk, times2saveidx(i)] = min(abs(EEG.times - times2save(i)));
end

pts_in_timewin = round((timewin/1000)*EEG.srate);

% create hann taper
hann_win = 0.5 * (1 - cos(2 * pi * (0:pts_in_timewin - 1)/(pts_in_timewin - 1)));

% freqs to save 
freqs2save = linspace(0, EEG.srate/2, floor(pts_in_timewin/2) + 1);

% result matrix
times2freqs = zeros(length(freqs2save), length(times2save));

for ti = 1:length(times2saveidx)
    start_idx = times2saveidx(ti) - floor(pts_in_timewin/2);
    % if pts_in_timewin is even, then needs to minus 1 caused by the central time.
    end_idx = times2saveidx(ti) + floor(pts_in_timewin/2) - mod(pts_in_timewin + 1, 2);
    
    windat = squeeze(EEG.data(channel2useidx, start_idx:end_idx, :));
    
    taperdat = bsxfun(@times, windat, hann_win');

    % fft(x, t, dim), make sure fft is over time
    fftdat = fft(taperdat, [], 1);

    times2freqs(:, ti) = mean(abs(fftdat(1:floor(pts_in_timewin/2) + 1, :)).^2, 2);
end 


%% Hilbert
filter_interval = 2;
filtered_dat = eegfilt(EEG.data(channel2useidx, :, :), EEG.srate, freq2plot - filter_interval, freq2plot + filter_interval);
hilbert_res = hilbert(reshape(filtered_dat, 640, 99));

%% Wavelet
t = 1;
times = linspace(-t, t, 640);
sin_wave = exp(2 * pi * 1i * freq2plot .* times);
std_gauss = 3/(2*pi*freq2plot);
exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
wavelet = sin_wave .* exp_wave;
wavelet_res = zeros(640, 99);
for tri = 1:99
    wavelet_res(:, tri) = conv(EEG.data(channel2useidx, :, tri), wavelet, 'same');
end
% wavelet_res = conv(mean(EEG.data(channel2useidx, :, :), 3), wavelet, 'same');

% plot

figure
hold on;
[junk, freq2plotidx] = min(abs(freqs2save - freq2plot));
plot(times2save, zscore(times2freqs(freq2plotidx, :)),'DisplayName','stFFT');
plot(times2save, zscore(mean(abs(hilbert_res(t_s:t_e, :).^2), 2)), 'DisplayName','filter-Hilbert');
% plot(times2save, zscore(abs(wavelet_res(t_s:t_e)).^2), 'DisplayName','wavelet');
plot(times2save, zscore(mean(abs(wavelet_res(t_s:t_e, :)).^2, 2)), 'g','DisplayName','wavelet');

xlabel('Time(ms)');
ylabel('Zscore power');
hold off
legend

