%{
2. Select one electrode from the scalp EEG dataset and convolve each wavelet with EEG data from all trials from that electrode. 

Apply the Matlab function real to the convolution result, as in convol 
(convol This will return the EEG data bandpass filtered at the peak frequency of the wavelet.
You learn more about why this is in the next chapter.
%}


frequencies = linspace(2, 30, 5);
times = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
eeg_data = EEG.data(47, :, 1);
num_img_freq = 3;

figure
for fi = 1:length(frequencies)
    sin_wave = exp(1i * 2 * pi * frequencies(fi) * times);
    exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
    morlet_wave = sin_wave .* exp_wave;
    
    subplot(5,num_img_freq,fi * num_img_freq - 2);
    plot(times, real(sin_wave));
    subplot(5,num_img_freq,fi * num_img_freq - 1);
    plot(times, morlet_wave);
    subplot(5,num_img_freq,fi * num_img_freq);
    hold on;
    plot(EEG.times, normalize(eeg_data, 'range'));
    plot(EEG.times, normalize(real(conv(eeg_data, morlet_wave, 'same')), 'range'));
    hold off;
end
