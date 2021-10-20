%{
1. Select one seed electrode and one frequency band and compute phase-based connectivity 
between that seed electrode and every other electrode. Use two methods for phase-based 
connectivity that were presented in this chapter, one that is volume conduction independent (e.g., PLI) 
and one that could produce spurious connectivity due to volume conduction (e.g., ISPC). 
Do not apply a baseline subtraction. 
Make topographical plots of seeded connectivity in a time window of your choice 
(e.g., 300-350 ms). 
What are the similarities and differences between results from the two methods, 
and what might be the reasons for the similarities and differences?

2. Now apply a baseline subtraction to the results (you can choose the baseline time period). 
Are there any changes in the plots after baseline subtraction 
(note that the color scaling will be different after baseline subtraction), 
and how do results from the two analyses compare with each other after baseline subtraction?

3. From the results in exercise 1 above, pick one "target" electrode 
(any electrode other than the seed) and provide evidence, using additional data analyses 
if necessary, for or against that measure of phase-based connectivity 
between that electrode and the seed being driven by volume conduction.

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

seed_channel = 'p1';
seed_freq = 5; % Hz
time2use = 325; 
timewin = 25;

baselinetm = [-400 -200];
time2useidx = dsearchn(EEG.times', time2use);
baselineidx   = dsearchn(EEG.times',baselinetm');
timewinidx = round(timewin*EEG.srate/1000);

time = -1:1/EEG.srate:1; 
s = 3.5/(2*pi*seed_freq);
wavelet = exp(2*1i*pi*seed_freq.*time) .* exp(-time.^2./(2*(s^2)));
half_of_wavelet_size = (length(time)-1)/2;

n_wavelet = length(time);
n_data = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

seedidx = find(strcmpi(seed_channel,{EEG.chanlocs.labels}));

wavelet_fft = fft(wavelet,n_convolution);
seed_fft = fft(reshape(EEG.data(seedidx, :, :), 1, n_data), n_convolution);
seed_conv_result = ifft(wavelet_fft.*seed_fft, n_convolution);
seed_conv_result = seed_conv_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
seed_sig = reshape(seed_conv_result, EEG.pnts, EEG.trials);

result = zeros(4, EEG.nbchan);

for chani=1:EEG.nbchan
    target_fft = fft(reshape(EEG.data(chani, :, :), 1, n_data), n_convolution);

    target_conv_result = ifft(wavelet_fft .* target_fft, n_convolution);
    target_conv_result = target_conv_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    target_sig = reshape(target_conv_result, EEG.pnts, EEG.trials);

    cdd = seed_sig .* conj(target_sig);

    seed_spec = sum(seed_sig .* conj(seed_sig), 2);
    target_spec = sum(target_sig .* conj(target_sig), 2);
    res_spec = sum(seed_sig .* conj(target_sig), 2);
    ic = abs(imag(res_spec./sqrt(seed_spec.*target_spec)));

    cdd_base = cdd(baselineidx(1):baselineidx(2), :);
    cdd = cdd(time2useidx-timewinidx:time2useidx+timewinidx, :);

    result(1, chani) = mean(abs(mean(exp(1i * angle(cdd)), 2)));
    
    result(2, chani) = mean(ic(time2useidx-timewinidx:time2useidx+timewinidx, :));

    result(3, chani) = mean(abs(mean(exp(1i * angle(cdd)), 2))) - mean(abs(mean(exp(1i * angle(cdd_base)), 2)));

    result(4, chani) = mean(ic(time2useidx-timewinidx:time2useidx+timewinidx, :)) - mean(ic(baselineidx(1):baselineidx(2), :));
end

subplot(2, 2, 1);
topoplot(result(1,:), EEG.chanlocs, 'style', 'both', 'maplimits', [0 0.8]);
title('ISPC-trials, no baseline correction');
cbar('vert',0);

subplot(2, 2, 2);
topoplot(result(2,:), EEG.chanlocs, 'style', 'both', 'maplimits', [0 0.15]);
title('Imaginary coherence, no baseline correction');
cbar('vert',0);

subplot(2, 2, 3);
topoplot(result(3,:), EEG.chanlocs, 'style', 'both', 'maplimits', [-0.1 0.1]);
title('ISPC-trials, baseline correction');
cbar('vert',0);

subplot(2, 2, 4);
topoplot(result(4,:), EEG.chanlocs, 'style', 'both', 'maplimits', [-0.1 0.1]);
title('Imaginary coherence, baseline correction');
cbar('vert',0);