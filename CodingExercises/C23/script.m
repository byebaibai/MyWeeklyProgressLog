%{

1. Perform PCA on broadband data using two time windows, 
one before and one after trial onset (e.g., - 500 to 0 ms and 100 to 600 ms).

2. Plot topographical maps and time courses of the first four components. 
To construct the PCA time courses, multiply the PCA weights defined by the pre - and posttrial time windows with the electrode 
time courses from the entire trial. Do you notice any differences in the topographical maps or time courses from before 
versus after stimulus onset? How would you interpret differences and/or similarities?

3. Repeat this exercise but after bandpass filtering in two different frequency bands. 
Make sure there are no edge artifacts in the pretrial time window (consider using reflection, 
if necessary, as described in figure 7.3). Justify your decision of frequency bands and time window width(s). 
    Comment on any qualitative similarities and differences you observe between frequency bands and time windows and 
    similarities and differences between the frequencyband-specific and broadband signal from the results obtained in 
    the previous exercise.

%}

% Same with Solution

load ../../data/sampleEEGdata.mat;

% exercise 1 - 2

pre_stim_time = [-500 0];
[~, preidx(1)] = min(abs(EEG.times - pre_stim_time(1)));
[~, preidx(2)] = min(abs(EEG.times - pre_stim_time(2)));

post_stim_time = [100 600];
[~, postidx(1)] = min(abs(EEG.times - post_stim_time(1)));
[~, postidx(2)] = min(abs(EEG.times - post_stim_time(2)));


% filter data
center_freq = 15; % in Hz
filter_frequency_spread  = 3; % Hz +/- the center frequency
transition_width  = 0.2;

% construct filter kernel
nyquist       = EEG.srate/2;
filter_order  = round(3*(EEG.srate/(center_freq-filter_frequency_spread)));

ffrequencies  = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread) (center_freq-filter_frequency_spread) (center_freq+filter_frequency_spread) (1+transition_width)*(center_freq+filter_frequency_spread) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

filter_result = filtfilt(filterweights,1,double(reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials))')';
filter_result = reshape(filter_result,EEG.nbchan,EEG.pnts,EEG.trials);


figure

covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = squeeze(filter_result(:, preidx(1):preidx(2), i)) - repmat(squeeze(mean(filter_result(:,preidx(1):preidx(2),i), 2)), 1, preidx(2) - preidx(1) + 1);
    covar = covar + (eeg*eeg')./(preidx(2) - preidx(1));
end
covar = covar ./ EEG.trials;
[pc_pre, eigvals] = eig(covar);
pc_pre = pc_pre(:, end:-1:1);

ylims = [[-10 10];[-4 4]; [-2 2]; [-3 3]];
for i=1:4
    subplot(2, 2, i);
    plot(EEG.times, pc_pre(:, i)' * squeeze(mean(filter_result, 3)));
    hold on;
    set(gca, 'xlim', [-300 1000], 'ylim', ylims(i, :));
    title([ 'PC #' num2str(i) ' time course']);
end

covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = squeeze(filter_result(:, postidx(1):postidx(2), i)) - repmat(squeeze(mean(filter_result(:,postidx(1):postidx(2),i), 2)), 1, postidx(2) - postidx(1) + 1);
    covar = covar + (eeg*eeg')./(postidx(2) - postidx(1));
end
covar = covar ./ EEG.trials;
[pc_post, eigvals] = eig(covar);
pc_post = pc_post(:, end:-1:1);

for i=1:4
    subplot(2, 2, i);
    plot(EEG.times, pc_post(:, i)' * squeeze(mean(filter_result, 3)));
    hold on;
    set(gca, 'xlim', [-300 1000], 'ylim', ylims(i, :));
    title([ 'PC #' num2str(i) ' time course']);
    legend({'pre-stim';'post-stim'})
end

figure

for i=1:4
    subplot(4, 2, i*2 - 1);
    topoplot(pc_pre(:, i), EEG.chanlocs, 'electrodes', 'off', 'plotrad', 0.53);
    title([ 'pre-stim PC #' num2str(i) ' map']);
    subplot(4, 2, i*2);
    topoplot(pc_post(:, i), EEG.chanlocs, 'electrodes', 'off', 'plotrad', 0.53);
    title([ 'post-stim PC #' num2str(i) ' map']);
    
end
