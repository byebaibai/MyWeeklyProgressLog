%{
1. Create a family of complex Morlet wavelets, ranging in frequencies from 2 Hz to 30 Hz in five steps.
%}
frequencies = linspace(2, 30, 5);
wavelets = zeros(5, 2001);
ts = linspace(-1, 1, 2001);
std_gauss = 5/(2*pi*30);

figure
for fi = 1:length(frequencies)
    sin_wave = exp(2 * pi * i * frequencies(fi) .* ts);
    exp_wave = exp(- ts.^2 / (2 * std_gauss.^2));
    wavelets(fi, :) = sin_wave .* exp_wave;
    subplot(5,1,fi);
    plot(ts, wavelets(fi, :));
end