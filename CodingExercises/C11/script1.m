%{
Reproduce the top three panels of figure 11.12 three times. 
First, perform time-domain convolution using the "manual" convolution method shown in chapter 10. 
Second, perform frequency-domain convolution using the discrete time Fourier transform 
    presented at the beginning of this chapter. 
Finally, perform frequency-domain convolution using the Matlab functions fftandifft (do not use the function cony). 
(You can optionally reproduce the bottom panel of figure 11.12 for the frequency domain analyses; 
keep in mind that the poualwer scaling is for display purposes only.)
%}

load ../data/sampleEEGdata.mat;

% time-domain convolution using the "manual" convolution method
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

figure
subplot(3, 1, 1);
plot(EEG.times, eeg_data);

subplot(3, 1, 2);
plot(EEG.times, gauss);

subplot(3, 1, 3);
hold on
plot(EEG.times, eeg_data);

padding_eeg_data = [zeros(1, length(gauss) - 1) eeg_data zeros(1, length(gauss) - 1)];
half_kernel_size = ceil((length(gauss) - 1)/2);
conv_result = zeros(1, length(t) + length(gauss) - 1);
for index = 1:length(conv_result) - half_kernel_size
    conv_result(index) = sum(gauss(end:-1:1) .* padding_eeg_data(index:index + length(gauss) - 1));
end
conv_result = conv_result(half_kernel_size: end - half_kernel_size);
plot(EEG.times, conv_result);

hold off

% frequency-domain convolution using the discrete time Fourier transform
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

figure
subplot(3, 1, 1);
plot(EEG.times, eeg_data);

subplot(3, 1, 2);
plot(EEG.times, gauss);

subplot(3, 1, 3);
hold on
plot(EEG.times, eeg_data);

half_kernel_size = ceil((length(gauss) - 1)/2);
padding_signal=[eeg_data zeros(1, length(gauss) - 1)];
padding_kernel=[gauss zeros(1, length(gauss) - 1)];
times = (0:length(padding_signal) - 1)/length(padding_signal);

eeg_fourier = zeros(1, length(padding_signal));
for index = 1:length(eeg_fourier)
    sin_wave = exp( -1i * 2 * pi * (index - 1) .* times);
    eeg_fourier(index) = sum(padding_signal .* sin_wave);
end
gauss_fourier = zeros(1, length(padding_kernel));
for index = 1:length(gauss_fourier)
    sin_wave = exp( -1i * 2 * pi * (index - 1) .* times);
    gauss_fourier(index) = sum(padding_kernel .* sin_wave);
end

result = ifft(eeg_fourier .* gauss_fourier);
result = result(half_kernel_size:end-half_kernel_size);
plot(EEG.times, result);
plot(EEG.times, result);
hold off;

% perform frequency-domain convolution using the Matlab functions fft and ifft
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

figure
subplot(3, 1, 1);
plot(EEG.times, eeg_data);

subplot(3, 1, 2);
plot(EEG.times, gauss);

subplot(3, 1, 3);
hold on
plot(EEG.times, eeg_data);

half_kernel_size = ceil((length(gauss) - 1)/2);
padding_signal=[eeg_data zeros(1, length(gauss) - 1)];
padding_kernel=[gauss zeros(1, length(gauss) - 1)];
result = ifft(fft(padding_signal) .* fft(padding_kernel));
result = result(half_kernel_size:end-half_kernel_size);
plot(EEG.times, result);

% plot(EEG.times, eeg_fourier);

hold off;