clc; close all; clear all;

%% Load EEG File

edf_file = '/Users/joeyroberts/Desktop/CAPSTONE/EEG2_Trim.edf'; % EEG file path
[hdr, record] = edfread(edf_file);
duration_minutes = height(hdr) / 60;

% Extract signal labels
full_signal_labels = hdr.Properties.VariableNames;
num_channels = length(full_signal_labels);

%% Prompt User to Select a Channel for Analysis

fprintf('Available EEG Channels:\n');
for i = 2:num_channels
    fprintf('%d: %s\n', i, full_signal_labels{i});
end

channel_idx = input('Enter the number of the EEG channel to analyze: ');
if channel_idx < 2 || channel_idx > num_channels
    error('Invalid channel selection. Please choose a valid channel number from the list.');
end

signal_data = hdr{:, channel_idx};
fprintf('Analyzing channel: %s\n', full_signal_labels{channel_idx});

ref_data = [hdr{:, 2:23}]; 

%% Estimate Sample Rate and Create Time Vector

sample_rate = 512;
time_limit = duration_minutes * 60;
sample_limit = round(time_limit * sample_rate);
time_per_second = (0:511) / sample_rate;

%% Align Signal Data with Time

ref_data_unpack = cell2mat(ref_data);
aligned_ref_data = reshape(ref_data_unpack, [], 22);

data_unpack = cell2mat(signal_data);
aligned_signal_data = reshape(data_unpack, [], 1);

aligned_time_data = linspace(600, time_limit+600, sample_limit);

% filtered_signal = aligned_signal_data;

%% Re-referencing

ref = mean(aligned_ref_data, 2);
re_ref_signal = aligned_signal_data - ref;

%% Apply Band-Pass Filtering (0.1 - 40 Hz)

low_cutoff = 0.1;  
high_cutoff = 40;  
order = 4;  % Filter order
nyquist = sample_rate / 2;  
wn = [low_cutoff, high_cutoff] / nyquist;  % Normalized cutoff frequencies

% Design a band-pass Butterworth filter
[b, a] = butter(order, wn, 'bandpass');

% Apply zero-phase filtering
filtered_signal = filtfilt(b, a, re_ref_signal);

%% Subplot Filtered Signal in 20-Second Intervals

seconds_to_test = 20;
num_intervals = 6;
starting_sec = 1370;
starting_index = (starting_sec-600)*sample_rate;
start_times = linspace(starting_index, starting_index + (num_intervals - 1) * seconds_to_test * sample_rate, num_intervals);

figure;
for i = 1:num_intervals
    start = round(start_times(i));
    finish = start + seconds_to_test * sample_rate;
    subplot(2, 3, i);
    plot(aligned_time_data(start:finish), filtered_signal(start:finish), 'r');
    ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
    title(['Filtered EEG Signal at ', num2str(round((start)/512) + 600), ' Seconds']);
    grid on;
    xlabel('Time (s)');
end

%% Phase Space Reconstruction

tau = 512;
figure;
for i = 1:num_intervals
    start = round(start_times(i));
    finish = start + seconds_to_test * sample_rate;
    X = filtered_signal(start:finish)';
    X3 = X(1:end - 2*tau);
    X2 = X(1+tau:end - tau);
    X1 = X(1+2*tau:end);
    subplot(2, 3, i);
    plot3(X1, X2, X3, 'b');
    xlim([-200 400]);
    ylim([-200 400]);
    zlim([-200 400]);
    grid on;
    xlabel('X(t)'); ylabel(['X(t-', num2str(tau/512), 'sec)']); zlabel(['X(t-', num2str(2*(tau/512)), 'sec)']);
    title(['Phase Space Reconstruction at ', num2str(round((start)/512) + 600), ' Seconds']);
end
