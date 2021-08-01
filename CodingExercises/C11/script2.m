%{
From the three sets of Matlab code you have for reproducing figure 11.12, 
run a computation time test. 
That is, time how long it takes Matlab to perform 1000 repetitions of each of the three methods 
for computing convolution that you generated in the previous exercise (do not plot the results each time). 
You can use the Matlab function pairs tic and toc to time a Matlab process. Plot the results in a bar plot, similar to figure 11.8.
%}

load ../data/sampleEEGdata.mat;

% time-domain convolution using the "manual" convolution method
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

tic;
for epoch = 1:1000
    padding_eeg_data = [zeros(1, length(gauss) - 1) eeg_data zeros(1, length(gauss) - 1)];
    half_kernel_size = ceil((length(gauss) - 1)/2);
    conv_result = zeros(1, length(t) + length(gauss) - 1);
    for index = 1:length(conv_result) - half_kernel_size
        conv_result(index) = sum(gauss(end:-1:1) .* padding_eeg_data(index:index + length(gauss) - 1));
    end
    conv_result = conv_result(half_kernel_size: end - half_kernel_size);
end
manual_conv_time = toc;
% frequency-domain convolution using the discrete time Fourier transform
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

tic;
for epoch = 1:1000
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
end
manual_fft_time = toc;

% perform frequency-domain convolution using the Matlab functions fft and ifft
t = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);
gauss = exp(- t.^2 / (2 * std_gauss^2))/30;
eeg_data = EEG.data(47, :, 1);

tic;
for epoch = 1:1000
    half_kernel_size = ceil((length(gauss) - 1)/2);
    padding_signal=[eeg_data zeros(1, length(gauss) - 1)];
    padding_kernel=[gauss zeros(1, length(gauss) - 1)];
    result = ifft(fft(padding_signal) .* fft(padding_kernel));
    result = result(half_kernel_size:end-half_kernel_size);
end
matlab_fft_time = toc;

cates = categorical({'Manual Conv','Manual FFT','Matlab FFT'});
cates = reordercats(cates, {'Manual Conv','Manual FFT','Matlab FFT'});
bar(cates, [manual_conv_time, manual_fft_time, matlab_fft_time])