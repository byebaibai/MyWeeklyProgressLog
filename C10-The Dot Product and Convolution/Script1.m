%{
    Problem Intro.
    Create two kernels for convolution: one that looks like an inverted U and one that looks like a decay function. 
    There is no need to be too sophisticated in generating, for example, a Gaussian and an exponential; 
    numerical approximations are fine
%}

kt = linspace(0, 4, 5)
kernel_decay = exp(-kt)
kernel_inverted_U = gaussmf(kt, [2 2])

hold on

plot(kernel_decay, 'DisplayName', "decay kernel");
plot(kernel_inverted_U, 'DisplayName', "inverted U kernel");
set(get(gca, 'XLabel'), 'String', 'Time');
set(get(gca, 'YLabel'), 'String', 'Voltage');
set(get(gca, 'Title'), 'String', 'Two Kernels for Convolution');

legend

hold off
