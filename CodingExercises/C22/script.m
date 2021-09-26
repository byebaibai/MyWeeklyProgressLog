%{

1. Based on the topographical maps of ERPs in figure 22.6 (plate 12), 
select one electrode whose activity you think might look similar before and after computing the surface Laplacian, 
and one electrode whose activity you think might look different before and after the surface Laplacian.

2. Perform a time-frequency decomposition of the data from those two electrodes both before and after computing the surface Laplacian 
(that is, compute the surface Laplacian on the raw data before applying a time-frequency decomposition). 
Compute both power (decibels normalized using a baseline period of your choice) and ITPC.

3. Plot the results using the same color scaling for before and after the surface Laplacian.

4. Are there any salient differences in the time-frequency power or ITPC results 
before versus after application of the surface Laplacian, and do the differences depend on the frequency? 
How would you interpret similarities and differences at different frequency bands?

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

% script 1
chans2plot = ["fc1" "po8"];
chans2plotidx = zeros(length(chans2plot), 1);
for chi = 1:length(chans2plot)
    chans2plotidx(chi) = find(strcmpi(chans2plot(chi), {EEG.chanlocs.labels}));
end

% script 2
% extract XYZ coordinates from EEG structure
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];


ori_data = EEG.data(chans2plotidx, :, :);
lap_data = laplacian_perrinX(EEG.data, X, Y, Z);
lap_data = lap_data(chans2plotidx, :, :);

window_time = [-300 1000];
[~, windowidx(1)] = min(abs(EEG.times - window_time(1)));
[~, windowidx(2)] = min(abs(EEG.times - window_time(2)));

baseline_time = [-400 -100];
[~, baseidx(1)] = min(abs(EEG.times - baseline_time(1)));
[~, baseidx(2)] = min(abs(EEG.times - baseline_time(2)));

freqs2plot = 2:1:40;

freq_interval = 1;

ori2filt = reshape(ori_data, size(ori_data, 1), []);
lap2filt = reshape(lap_data, size(lap_data, 1), []);

[ori_db, ori_itpc] = getHilbertResult(ori_data, chans2plot, freqs2plot, baseidx, freq_interval, EEG.srate);
[lap_db, lap_itpc] = getHilbertResult(lap_data, chans2plot, freqs2plot, baseidx, freq_interval, EEG.srate);

for chi = 1:length(chans2plot)
    figure
    subplot(2, 2, 1)
    contourf(EEG.times(windowidx(1):windowidx(2)), freqs2plot, squeeze(ori_db(chi, :, windowidx(1):windowidx(2))), 40, 'linecolor','none')
    colormap('jet')
    set(gca, 'YDir', 'normal', 'clim', [-3 3],'xlim',[-300 1000]);
    title(chans2plot(chi) + ' volt. power');
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');

    subplot(2, 2, 2)
    contourf(EEG.times(windowidx(1):windowidx(2)), freqs2plot, squeeze(ori_itpc(chi, :, windowidx(1):windowidx(2))), 40, 'linecolor','none')
    colormap('jet')
    set(gca, 'YDir', 'normal', 'clim', [.05 .25],'xlim',[-300 1000]);
    title(chans2plot(chi) + ' volt. itpc');
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');

    subplot(2, 2, 3)
    contourf(EEG.times(windowidx(1):windowidx(2)), freqs2plot, squeeze(lap_db(chi, :, windowidx(1):windowidx(2))), 40, 'linecolor','none')
    colormap('jet')
    set(gca, 'YDir', 'normal', 'clim', [-3 3],'xlim',[-300 1000]);
    title(chans2plot(chi) + ' Lap. power');
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');

    subplot(2, 2, 4)
    contourf(EEG.times(windowidx(1):windowidx(2)), freqs2plot, squeeze(lap_itpc(chi, :, windowidx(1):windowidx(2))), 40, 'linecolor','none')
    colormap('jet')
    set(gca, 'YDir', 'normal', 'clim', [.05 .25],'xlim',[-300 1000]);
    title(chans2plot(chi) + ' Lap. itpc');
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
end

function [db, itpc]=getHilbertResult(data, chans2plot, freqs2plot, baseidx, freq_interval, srate)
    data2filt = reshape(data, size(data, 1), []);
    hilbert_res = zeros(length(chans2plot), length(freqs2plot), size(data, 2),  size(data, 3));
    for fi = 1:length(freqs2plot)
        freq2plot = freqs2plot(fi);
        filtered_dat = reshape(eegfilt(data2filt, srate, freq2plot-freq_interval, freq2plot+freq_interval), size(data));
        for chi = 1:length(chans2plot)
            hilbert_res(chi, fi, :, :) = hilbert(squeeze(filtered_dat(chi, :, :)));
        end
    end
    hilbert_power = squeeze(mean(abs(hilbert_res).^2, 4));
    hilbert_baseline_power = squeeze(mean(hilbert_power(:, :, baseidx(1):baseidx(2)), 3));
    db = 10 * log10(bsxfun(@rdivide, hilbert_power, hilbert_baseline_power));
    itpc = squeeze(abs(mean(exp(1i * angle(hilbert_res)), 4)));
end