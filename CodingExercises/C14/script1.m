%{

Pick two frequencies (e.g., 5 Hz and 25 Hz) and one electrode and perform complex Morlet wavelet convolution 
and filter-Hilbert using those two frequencies as the peak/ center frequencies for all trials. 
Plot the resulting power and the bandpass-filtered signal (that is, the real component of the analytic signal) 
from each method. Plot one single trial (you can choose the trial randomly but plot the same trial for both methods) 
and then plot the average of all trials. Describe some similarities and differences between the results of 
the two time-frequency decomposition methods.

%}

load ../../data/sampleEEGdata.mat;
eegdata = EEG.data;
times = linspace(-3, 3, 640);
selected_electrode = 47;
selected_trial = 1;
frequencies = [2 25];
wavelets = zeros(length(frequencies), 640);
wavelet_frequency_spread = 3; % number of wavelet cycles


eeg_from_one_trial = squeeze(eegdata(selected_electrode, :, selected_trial));
eeg_from_all_trial = squeeze(mean(eegdata(selected_electrode, :, :), 3));

% s =  3./(2*pi*frex);
% s = 10./(2*pi*frex);

t_s = find(EEG.times >= -500, 1);
t_e = find(EEG.times <= 1000, 1, 'last');
figure
for fi = 1:length(frequencies)
    
    sin_wave = exp(2 * pi * 1i * frequencies(fi) .* times);
    gau_wave = exp(-times.^2./(2*(wavelet_frequency_spread/(2*pi*frequencies(fi)))^2));

    wavelets(fi, :) = sin_wave .* gau_wave;

    filtered_result = eegfilt(eeg_from_all_trial, EEG.srate, frequencies(fi) - 2, frequencies(fi) + 2);
    wavelet_result = conv(filtered_result, wavelets(fi, :), 'same');
    hilbert_result = hilbert(filtered_result')'; % time should be in the first dimension. 
    
    subplot(2, 1, fi)
    hold on;
    plot(EEG.times(t_s:t_e), abs(wavelet_result(t_s:t_e)));
    plot(EEG.times(t_s:t_e), abs(hilbert_result(t_s:t_e)));
    hold off;
end