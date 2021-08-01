%{
    Problem Intro.
    Convolve these two kernels with 50 time points of EEG data from one electrode. 
    Make a plot showing the kernels, the EEG data, and the result of the convolution between the data and each kernel. 
    Use time-domain convolution as explained in this chapter and as illustrated in the online Matlab code. 
    Based on visual inspection, what is the effect of convolving the EEG data with these two kernels?
%}


load ../data/sampleEEGdata.mat;

ks = 5;
num_dot = 50;

t = linspace(0, num_dot-1, num_dot);
eeg_data = EEG.data(1, 1:num_dot, 1);

kt = linspace(0, ks-1, ks);
kernel_decay = exp(-kt);
kernel_inverted_U = gaussmf(kt, [2 2]);

result_decay = conv(eeg_data, kernel_decay);
result_inverted_U = conv(eeg_data, kernel_inverted_U);
result_scaled_decay = conv(eeg_data, kernel_decay)/sum(kernel_decay);
result_scaled_inverted_U = conv(eeg_data, kernel_inverted_U)/sum(kernel_inverted_U);

padding_eeg_data = [zeros(1, 4) eeg_data zeros(1, 4)];
result_padding_inverted_U = zeros(1, num_dot + ks - 1);
for index = 1:num_dot + ks - 1
    result_padding_inverted_U(index) = dot(padding_eeg_data(1, index:index+4), kernel_inverted_U);
end
result_scaled_padding_inverted_U = result_padding_inverted_U/sum(kernel_inverted_U);

hold on

plot(eeg_data, 'DisplayName', "EEG Data");
plot(result_scaled_decay, 'DisplayName', "Scaled Decay Result");
plot(result_scaled_inverted_U, 'DisplayName', "Scaled Inverted U Result", "LineWidth", 2);
plot(result_scaled_padding_inverted_U, 'DisplayName', "Manual Scaled Padding Inverted U Result");
set(get(gca, 'XLabel'), 'String', 'Time');
set(get(gca, 'YLabel'), 'String', 'Voltage');


legend;

hold off