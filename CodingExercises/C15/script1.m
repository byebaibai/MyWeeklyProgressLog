%{

Compute the short-time FFT at each electrode and make topographical maps of thetaband (around 6 Hz) power 
and alpha-band (around 10 Hz) power at 150 ms and 700 ms.

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

timewin = 400; % in ms, for short-time FFT
times2save = -300:50:1000; % in ms
channel2plot = 47;
frequencies2plot = [6 10]; % in Hz
timepoints2plot = [150 700]; % ms
channel2useidx = (linspace(1, 64, 64) == channel2plot);

times2saveidx = zeros(size(times2save));
for i = 1:length(times2save)
    [junk, times2saveidx(i)] = min(abs(EEG.times - times2save(i)));
end

pts_in_timewin = round((timewin/1000)*EEG.srate);

% create hann taper
hann_win = 0.5 * (1 - cos(2 * pi * (0:pts_in_timewin - 1)/(pts_in_timewin - 1)));
hann_win_all = repmat(hann_win, 64, 1);
% freqs to save 
freqs2save = linspace(0, EEG.srate/2, floor(pts_in_timewin/2) + 1);

% result matrix
times2freqs = zeros(64, length(freqs2save), length(times2save));

for ti = 1:length(times2saveidx)
    start_idx = times2saveidx(ti) - floor(pts_in_timewin/2);
    % if pts_in_timewin is even, then needs to minus 1 caused by the central time.
    end_idx = times2saveidx(ti) + floor(pts_in_timewin/2) - mod(pts_in_timewin + 1, 2);
    
    windat = squeeze(EEG.data(:, start_idx:end_idx, :));
    
    taperdat = bsxfun(@times, windat, hann_win_all);

    % fft(x, t, dim), make sure fft is over time
    fftdat = fft(taperdat, [], 2);

    times2freqs(:, :, ti) = mean(abs(fftdat(:, 1:floor(pts_in_timewin/2) + 1, :)).^2, 3);
end 

for fi = 1:length(frequencies2plot)
    [junk, freq2plotidx] = min(abs(freqs2save - frequencies2plot(fi)));
    for ti = 1:length(timepoints2plot)
        [junk, time2plotidx] = min(abs(times2save - timepoints2plot(ti)));
        subplot(2, 2, ti * 2 + fi - 2)
        titleStr = sprintf('%dms, %dHz',timepoints2plot(ti),frequencies2plot(fi));
        title(titleStr)
        topoplot(times2freqs(:, freq2plotidx, time2plotidx), EEG.chanlocs);
    end
end

