%% Figure 13.1

% parameters...
srate = 500; % sampling rate in Hz
f     = 10; % frequency of wavelet in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
s     = 6/(2*pi*f);

% and together they make a wavelet
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 

figure
subplot(221)
% show the projection onto the real axis
plot3(time,real(wavelet),imag(wavelet),'m')
xlabel('Time (ms)'), ylabel('real axis')
view(0,90)
title('Projection onto real and time axes')

% show the projection onto the imaginary axis
subplot(222)
plot3(time,real(wavelet),imag(wavelet),'g')
xlabel('Time (ms)'), ylabel('imaginary axis')
view(0,0)
title('Projection onto imaginary and time axes') 
 
% plot projection onto real and imaginary axes
subplot(223)
plot3(time,real(wavelet),imag(wavelet),'k')
ylabel('real axis'), zlabel('imag axis')
view(90,0)
title('Projection onto imaginary and time axes')

% plot real and imaginary projections simultaneously
subplot(224)
plot(time,real(wavelet),'b')
hold on
plot(time,imag(wavelet),'b:')
legend({'real part';'imaginary part'})