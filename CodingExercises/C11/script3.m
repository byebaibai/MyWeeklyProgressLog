%{
Generate a time series by creating and summing sine waves, as in figure 11.2B. [Done]
Use between two and four sine waves, so that the individual sine waves are still somewhat visible in the sum. [Done]
Perform a Fourier analysis (you can use the fft function) on the resulting time series and plot the power structure. 
Confirm that your code is correct by comparing the frequencies with nonzero power to the frequencies of the sine waves that you generated. 
Now try adding random noise to the signal before computing the Fourier transform. 
    First, add a small amount of noise so that the sine waves are still visually recognizable. 
    Next, add a large amount of noise so that the sine waves are no longer visually recognizable in the time domain data. 
    Perform a Fourier analysis on the two noisy signals and plot the results. 
What is the effect of a small and a large amount of noise in the power spectrum? 
Are the sine waves with noise easier to detect in the time domain or in the frequency domain, 
or is it equally easy/difficult to detect a sine wave in the presence of noise?
%}
 
ts = linspace(0, 1, 1001);
A = 10;
sin1 = A * sin(2 * pi * 40 * ts);
sin2 = A * sin(2 * pi * 20 * ts);
sin3 = A * sin(2 * pi * 10 * ts);
sin4 = A * sin(2 * pi * 5 * ts);
sin_wave = sin1 + sin2 + sin3 + sin4;
noise_small = 50 * rand(size(ts));
noise_large = 100 * rand(size(ts));
subplot(421);
plot(ts, sin1);

subplot(422);
plot(ts, sin2);

subplot(423);
plot(ts, sin3);

subplot(424);
plot(ts, sin4);

subplot(4, 2,[5, 6]);
plot(ts, sin_wave);

subplot(4, 2,[7, 8]);
fft_sin_wave = fft(sin_wave)/length(sin_wave);
bar(abs(fft_sin_wave(2:100)*2));

figure;
sin_wave_small_noise = sin_wave + noise_small;
sin_wave_large_noise = sin_wave + noise_large;
subplot(221);
plot(ts, sin_wave_small_noise)
subplot(222);
fft_sin_wave_small_noise = fft(sin_wave_small_noise)/length(sin_wave_small_noise);
bar(abs(fft_sin_wave_small_noise(2:100)*2));

subplot(223);
plot(ts, sin_wave_large_noise)
subplot(224);
fft_sin_wave_large_noise = fft(sin_wave_large_noise)/length(sin_wave_large_noise);
bar(abs(fft_sin_wave_large_noise(2:100)*2));