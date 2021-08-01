%{
1. Create a family of Morlet wavelets ranging in frequency from 2 Hz to 30 Hz in five steps. 

Apply the Matlab function real to the convolution result, as in convol 
(convol This will return the EEG data bandpass filtered at the peak frequency of the wavelet.
You learn more about why this is in the next chapter.
%}

frequencies = linspace(2, 30, 5);
times = linspace(-1, 1, 2001);
std_gauss = 5/(2*pi*30);

figure
for fi = 1:length(frequencies)
    sin_wave = exp(1i * 2 * pi * frequencies(fi) * times);
    exp_wave = exp(-(times).^2 / (2 * std_gauss.^2));
    morlet_wave = real(sin_wave) .* exp_wave;
    
    subplot(5,2,fi * 2 - 1);
    plot(times, real(sin_wave));
    subplot(5,2,fi * 2);
    plot(times, morlet_wave);
end