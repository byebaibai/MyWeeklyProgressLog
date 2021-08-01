%{
Plot the broadband ERP (without any convolution). 
Thus, you will have six subplots in one figure. 
How do the wavelet-convolved ERPs compare with the broadband ERP? 
Are there dynamics revealed by the wavelet-convolved ERPs that are not apparent in the broadband ERP, 
and are there dynamics in the broadband ERP that are not apparent in the waveletconvolved ERPs? 
Base your answer on qualitative visual inspection of the results; statistics or other quantitative comparisons are not necessary.

%}

frequencies = linspace(2, 30, 5);
times = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
eeg_data = mean(EEG.data(47, :, :),3);
num_img_freq = 6;

figure
for fi = 1:length(frequencies)
    sin_wave = exp(1i * 2 * pi * frequencies(fi) * times);
    exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
    morlet_wave = sin_wave .* exp_wave;
    
    subplot(5,num_img_freq,fi * num_img_freq - 5);
    plot(times, sin_wave);

    subplot(5,num_img_freq,fi * num_img_freq - 4);
    plot(times, morlet_wave);

    subplot(5,num_img_freq,fi * num_img_freq - 3);
    fft_kernel = fft(morlet_wave)/length(morlet_wave);
    bar(abs(fft_kernel(2:100)) * 2);


    subplot(5,num_img_freq,fi * num_img_freq - 2);
    hold on;
    plot(EEG.times, normalize(eeg_data, 'range'));
    conv_signal = conv(eeg_data, morlet_wave, 'same');
    plot(EEG.times, normalize(real(conv_signal), 'range'));
    hold off;

    subplot(5,num_img_freq,fi * num_img_freq - 1);
    fft_signal = fft(eeg_data);
    bar(abs(fft_signal(2:100)) * 2);

    subplot(5,num_img_freq,fi * num_img_freq);
    conv_signal = conv(eeg_data, morlet_wave, 'same');
    fft_conv_signal = fft(conv_signal);
    bar(abs(fft_conv_signal(2:100)) * 2);

end
