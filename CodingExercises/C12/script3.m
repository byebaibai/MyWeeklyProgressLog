%{
3. Average the result of convolution over all trials and plot an ERP corresponding to each wavelet frequency. Each frequency should be in its own subplot.

%}

frequencies = linspace(2, 30, 5);
times = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
eeg_data = mean(EEG.data(47, :, :),3);
num_img_freq = 2;

figure
for fi = 1:length(frequencies)
    sin_wave = exp(1i * 2 * pi * frequencies(fi) * times);
    exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
    morlet_wave = sin_wave .* exp_wave;
    
    subplot(5,num_img_freq,fi * num_img_freq - 1);
    plot(times, morlet_wave);
    subplot(5,num_img_freq,fi * num_img_freq);
    hold on;
    
    plot(EEG.times, normalize(eeg_data, 'range'));
    plot(EEG.times, normalize(real(conv(eeg_data, morlet_wave, 'same')), 'range'));
    hold off;
end
