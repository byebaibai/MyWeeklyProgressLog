%{
3. Extract power and phase from the result of complex wavelet convolution and store in a time x frequency x electrodes x power/phase matrix 
(thus, a 640 x 5 x 64 x 2 matrix).
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
dataFrom1stTrial = eegdata(:, :, 98);
convResults = zeros(length(frequencies), size(dataFrom1stTrial, 1), size(dataFrom1stTrial, 2));
powerAndPhaseResults = zeros(size(dataFrom1stTrial, 2), length(frequencies), size(dataFrom1stTrial, 1), 2);

for fi = 1:length(frequencies)
    for elcs = 1:size(dataFrom1stTrial, 1)
        convResults(fi, elcs, :) = conv(wavelets(fi, :)', dataFrom1stTrial(elcs, :), 'same');
        powerAndPhaseResults(:, fi, elcs, 1) = convResults(fi, elcs, :) .* conj(convResults(fi, elcs, :));
        powerAndPhaseResults(:, fi, elcs, 2) = atan(imag(convResults(fi, elcs, :)) ./ real(convResults(fi, elcs, :)));
    end
    subplot(5, 2, fi * 2 - 1);
    plot(ts, squeeze(powerAndPhaseResults(:, fi, 20, 1)));
    subplot(5, 2, fi * 2);
    plot(ts, squeeze(powerAndPhaseResults(:, fi, 20, 2)));
end