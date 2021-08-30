%{

For each of these three electrodes, 
compute wITPCz using reaction time as the trial-varying modulator. 
Perform this analysis for all time-frequency points to generate time-frequency 
maps of the relationship between phase and reaction time. 
Do the time-frequency maps of wITPCz look different from the time-frequency 
maps of ITPC? Do you see any striking patterns in the ITPCz results, 
and do the results differ across the different electrodes (don't worry 
about statistics, base your judgment on qualitative patterns)? How would you 
interpret the results if they were statistically significant?

%}


load ../../data/sampleEEGdata.mat;

channels2plot = ["fz" "p8" "oz"];
channels2plotidx = zeros(length(channels2plot), 1);
for chi = 1:length(channels2plot)
    channels2plotidx(chi) = find(strcmpi(channels2plot(chi),{EEG.chanlocs.labels}));
end

t_start = -300;
t_end = 1000;
[junk, t_start_idx] = min(abs(EEG.times - t_start));
[junk, t_end_idx] = min(abs(EEG.times - t_end));

baseline_start = -300;
baseline_end = -100;
[junk, baseline_start_idx] = min(abs(EEG.times - baseline_start));
[junk, baseline_end_idx] = min(abs(EEG.times - baseline_end));

% initialize matrix to store RTs
rts = zeros(size(EEG.epoch));
for ei=1:EEG.trials
    % find which event is time=0, and take the latency of the event there after.
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
    % use try-catch in case of no response
    try
        rts(ei) = EEG.epoch(ei).eventlatency{time0event+1};
    catch me;
        rts{ei} = NaN;
    end
end

freqs2plot = 2:1:30;
% wavelet
% times = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
% wavelet_res = zeros(length(channels2plot), length(freqs2plot), size(EEG.data, 2),size(EEG.data, 3));
% for fi = 1:length(freqs2plot)
%     freq2plot = freqs2plot(fi);
%     std_gauss = 4/(2*pi*freq2plot);
%     exp_wave = exp( (-(times).^2) ./ (2 * (std_gauss.^2)) );
%     sin_wave = exp(2 * pi * 1i * freq2plot .* times);
%     wavelet = sin_wave .* exp_wave / freq2plot;
%     for chi = 1:length(channels2plot)
%         for tri = 1:size(EEG.data, 3)
%             wavelet_res(chi, fi, :, tri) = conv(EEG.data(chi, :, tri), wavelet, 'same');
%         end
%     end
% end

% wavelet_power = squeeze(mean(abs(wavelet_res).^2, 4));
% wavelet_baseline_power = squeeze(mean(wavelet_power(:, :, baseline_start_idx:baseline_end_idx), 3));
% wavelet_db = 10 * log10(bsxfun(@rdivide, wavelet_power, wavelet_baseline_power));


% filter-hilbert
freq_interval = 1; %Hz
data2filt = reshape(EEG.data, size(EEG.data, 1), []);
hilbert_res = zeros(length(channels2plot), length(freqs2plot), EEG.pnts, EEG.trials);
hilbert_witpc = zeros(length(channels2plot), length(freqs2plot), EEG.pnts);
hilbert_witpc_z = zeros(length(channels2plot), length(freqs2plot), EEG.pnts);

for fi = 1:length(freqs2plot)
    freq2plot = freqs2plot(fi);
    filtered_dat = reshape(eegfilt(data2filt, EEG.srate, freq2plot - freq_interval, freq2plot + freq_interval), size(EEG.data));
    for chi = 1:length(channels2plot)
        hilbert_res(chi, fi, :, :) = hilbert(squeeze(filtered_dat(channels2plotidx(chi, :), :, :)));

        hilbert_angle = squeeze(angle(hilbert_res(chi, fi, :, :)));
        hilbert_witpc(chi, fi, :) = squeeze(abs(mean(bsxfun(@times, rts, exp(1i * hilbert_angle)), 2)));
        % permutation testing
        perm_witpc = zeros(EEG.pnts, 1000);
        for i=1:1000
            perm_witpc(:, i) = abs(mean(rts(randperm(EEG.trials)) .* exp(1i * hilbert_angle), 2));
        end
        hilbert_witpc_z(chi, fi, :) = (squeeze(hilbert_witpc(chi, fi, :))-mean(perm_witpc, 2))./std(perm_witpc, 0, 2);
    end
end
hilbert_power = squeeze(mean(abs(hilbert_res).^2, 4));
hilbert_baseline_power = squeeze(mean(hilbert_power(:, :, baseline_start_idx:baseline_end_idx), 3));
hilbert_db = 10 * log10(bsxfun(@rdivide, hilbert_power, hilbert_baseline_power));
hilbert_itpc = squeeze(abs(mean(exp(1i * angle(hilbert_res)), 4)));


figure
for chi = 1:length(channels2plot)
    subplot(2, 3, chi);
    contourf(EEG.times(t_start_idx:t_end_idx), freqs2plot, squeeze(hilbert_witpc(chi, :, t_start_idx:t_end_idx)), 40, 'linecolor','none')
    colormap('jet')
    colorbar;
    set(gca, 'YDir', 'normal','xlim',[-300 1000]);
    title('wITPC from ' + channels2plot(chi));
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');

    subplot(2, 3, chi + 3);  
    contourf(EEG.times(t_start_idx:t_end_idx), freqs2plot, squeeze(hilbert_witpc_z(chi, :, t_start_idx:t_end_idx)), 40, 'linecolor','none')
    colormap('jet')
    colorbar;
    set(gca, 'YDir', 'normal','xlim',[-300 1000]);
    title('wITPCz from ' + channels2plot(chi));
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
end