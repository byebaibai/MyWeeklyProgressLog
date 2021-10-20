%{

1. Perform a power correlation analysis over time. 
Pick two electrodes and use a sliding time segment of three cycles 
(1.5 cycles on either side of each center time point). Average the results over trials. 
Perform this analysis at three frequencies and plot the time series of correlation coefficients. 
Next, repeat the analysis twice, using fixed time-segment lengths of 150 ms and 900 ms. 

Do the results differ according to the time segment length and the frequency band, 
and how are they different? In what situations would it be beneficial to use each window width parameter, 
and in what situations might problems or limitations arise? 

2. Select two "seed" time-frequency-electrode windows and perform an exploratory power 
correlation over trials at one selected "target" electrode, as in figure 27.6C (plate 19). 
Show the results in separate plots, and then show a time-frequency plot of correlation 
coefficient differences between the two seeds (Fisher-Z transform the coefficients before 
subtraction). Are there any striking qualitative differences between the two plots, 
and did plotting the difference map make the differences easier or more difficult to interpret?

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

% exercise 1

channel1 = 'fz';
channel2 = 'o1';

freqs2use = [5, 10, 15];

% wavelet parameters
time = -1:1/EEG.srate:1;
half_of_wavelet = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 2.3;

% FFT of data
fft_data1 = fft(reshape(EEG.data(strcmpi(channel1,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data2 = fft(reshape(EEG.data(strcmpi(channel2,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);

colors = ['b';'g';'r'];
for fi=1:length(freqs2use)

    freq = freqs2use(fi);
    three_times_win = round(1000/freq);
    windows = [150 900 three_times_win];
    % create wavelet and get its FFT
    fft_wavelet = fft(exp(2*1i*pi*freq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*freq))^2))/freq,n_convolution);
    
    % convolution for all three sites (save only power)
    convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet+1:end-half_of_wavelet);
    conv_result1 = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet+1:end-half_of_wavelet);
    conv_result2 = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;

    % downsample and rank transform all data
    conv1 = tiedrank(conv_result1')';
    conv2 = tiedrank(conv_result2')';

    subplot(3, 1, fi);
    hold on;
    for wi=1:length(windows)
        times2save = -300:20:800;
        times2saveidx = dsearchn(EEG.times',times2save');
        timewinidx = round(windows(wi)/2*EEG.srate/1000);
        corr_res = zeros(1, length(times2save));

        for ti=1:length(times2save)
            % compute bivariate correlations
            t = times2saveidx(ti);
            [junk, timeidx(1)] = min(abs(EEG.times-times2save(ti) + windows(wi)/2)); 
            [junk, timeidx(2)] = min(abs(EEG.times-times2save(ti) - windows(wi)/2));
            timeidx
            % r_st = 1 - 6*sum((conv1(t-timewinidx:t+timewinidx,:)-conv2(t-timewinidx:t+timewinidx,:)).^2, 2)/(EEG.trials*(EEG.trials^2-1));
            r_st = 1 - 6*sum((conv1(timeidx(1):timeidx(2),:)-conv2(timeidx(1):timeidx(2),:)).^2, 2)/(EEG.trials*(EEG.trials^2-1));
            
            corr_res(ti) = mean(r_st);
        end
        plot(times2save, corr_res, colors(wi));
        if wi==3 
            legend({[num2str(windows(1)) ' ms'];[num2str(windows(2)) ' ms'];['3f (' num2str(windows(3)) ' ms)']});
        end
    end
    title(['Correlation between ' channel1 ' and ' channel2 ' at ' num2str(freq) ' Hz']);
    set(gca,'xlim',[-300 800],'ylim',[-.2 .2])
    xlabel('Time (ms)');
    ylabel('Correlation coef.');
    hold off;
end
    