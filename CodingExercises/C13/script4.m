%{
4. Make topographical plots of power and phase at 180 ms at all frequencies 
(hint: you may need to use the squeeze function to remove singleton dimensions). 
Arrange the plots in one figure with five columns for frequency and two rows for power/phase. 
Put labels in the plot so it is clear which topographical maps correspond to which frequencies.
%}
frequencies = linspace(2, 30, 5);
wavelets = zeros(5, 640);
ts = linspace(-1, 1, 640);
std_gauss = 5/(2*pi*30);

for fi = 1:length(frequencies)
    sin_wave = exp(2 * pi * 1i * frequencies(fi) .* ts);
    exp_wave = exp(- ts.^2 / (2 * std_gauss.^2));
    wavelets(fi, :) = sin_wave .* exp_wave;
end

figure 
load ../../data/sampleEEGdata.mat;
eegdata = EEG.data;
dataFrom1stTrial = eegdata(:, :, 98);
convResults = zeros(length(frequencies), size(dataFrom1stTrial, 1), size(dataFrom1stTrial, 2));
powerAndPhaseResults = zeros(size(dataFrom1stTrial, 2), length(frequencies), size(dataFrom1stTrial, 1), 2);

time = find((179 < EEG.times)&(EEG.times < 181), 1);
for fi = 1:length(frequencies)
    for elcs = 1:size(dataFrom1stTrial, 1)
        convResults(fi, elcs, :) = conv(wavelets(fi, :), dataFrom1stTrial(elcs, :), 'same');
        powerAndPhaseResults(:, fi, elcs, 2) = atan(imag(convResults(fi, elcs, :)) ./ real(convResults(fi, elcs, :)));
        powerAndPhaseResults(:, fi, elcs, 1) = convResults(fi, elcs, :) .* conj(convResults(fi, elcs, :));
    end

    subplot(2, 5, fi);
    topoplot(squeeze(powerAndPhaseResults(time, fi, :, 1)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');
    subplot(2, 5, fi + 5);
    topoplot(squeeze(powerAndPhaseResults(time, fi, :, 2)), EEG.chanlocs, 'style', 'map', 'electrodes', 'off');
end