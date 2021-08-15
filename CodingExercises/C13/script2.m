%{
2. Convolve each wavelet with EEG data from all electrodes and from one trial.
%}
frequencies = linspace(2, 30, 5);
wavelets = zeros(5, 640);
ts = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);

for fi = 1:length(frequencies)
    sin_wave = exp(2 * pi * i * frequencies(fi) .* ts);
    exp_wave = exp(- ts.^2 / (2 * std_gauss.^2));
    wavelets(fi, :) = sin_wave .* exp_wave;
end

figure 
load ../../data/sampleEEGdata.mat;
eegdata = EEG.data;
dataFrom1stTrial = eegdata(:, :, 1);
convResults = zeros(length(frequencies), size(dataFrom1stTrial, 1), size(dataFrom1stTrial, 2));

for fi = 1:length(frequencies)
    for elcs = 1:size(dataFrom1stTrial, 1)
        convResults(fi, elcs, :) = conv(wavelets(fi, :)', dataFrom1stTrial(elcs, :), 'same');
    end
    subplot(5, 2, fi * 2 - 1);
    plot(ts, real(squeeze(convResults(fi, 20, :))));
    subplot(5, 2, fi * 2);
    plot(ts, imag(squeeze(convResults(fi, 20, :))));
end
