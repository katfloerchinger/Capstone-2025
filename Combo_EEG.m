clc; close all; clear all;

%% Load EEG File

edf_file = 'C:\Users\sfmdf\OneDrive\Documents\MATLAB\Capstone\EEG2.edf'; % EDF file path
[hdr, record] = edfread(edf_file);
duration_minutes = height(hdr)/60;

% Extract signal labels (keeping all, but analysis starts at column 2)
full_signal_labels = hdr.Properties.VariableNames;  
num_channels = length(full_signal_labels);  % Total number of channels

%% Prompt User to Select a Channel for Analysis

fprintf('Available EEG Channels:\n');
for i = 2:num_channels  % Display all channels except the first one
    fprintf('%d: %s\n', i, full_signal_labels{i});
end

% Ask for user input
channel_idx = input('Enter the number of the EEG channel to analyze: ');

% Validate input
if channel_idx < 2 || channel_idx > num_channels
    error('Invalid channel selection. Please choose a valid channel number from the list.');
end

% Extract selected channel data
signal_data = hdr{:, channel_idx};
fprintf('Analyzing channel: %s\n', full_signal_labels{channel_idx});

%% Estimate Sample Rate and Create Time Vector

sample_rate = 512; 
time_limit = duration_minutes * 60; 
sample_limit = round(time_limit * sample_rate);  % Sample limit for 2 minutes

% Create a time vector for plotting
time_per_second = (0:511) / sample_rate;  % Time vector for 512 points per second

%% Align Signal Data with Time (Before Omitting Any Data)

aligned_signal_data = [];
aligned_time_data = [];

for j = 1:time_limit
    % Extract the 512 data points per second
    second_data = signal_data{j};  
    second_time = time_per_second + (j-1);  

    % Append to the aligned signal data and time
    aligned_signal_data = [aligned_signal_data; second_data];
    aligned_time_data = [aligned_time_data; second_time'];
end

%% First Display the Raw EEG Signal Before Asking for Omission

figure;
plot(aligned_time_data, aligned_signal_data);
ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
title(['EEG Signal - ', full_signal_labels{channel_idx}]);
grid on;
xlabel('Time (s)');
fprintf('Look at the EEG plot above to determine how much time to omit.\n');

%% Ask How Much Startup Data to Omit

omit_time = input('Enter how many seconds of data to omit from the beginning: ');

% Convert time to number of data points
omit_samples = omit_time * sample_rate;

% Ensure valid input
if omit_samples < 0 || omit_samples >= length(aligned_signal_data)
    error('Invalid omission time. Must be within range of available data.');
end

% Remove first 'omit_samples' data points
aligned_signal_data(1:omit_samples) = [];
aligned_time_data(1:omit_samples) = [];

%% Identify and Compute Reference Channels
reference_channels = {'Fp1', 'Fp2', 'F7', 'F8','O1', 'O2'};  % T3 and T4, T5 and T6 not here
reference_indices = find(ismember(full_signal_labels, reference_channels));

% Ensure all reference channels are found
if length(reference_indices) < length(reference_channels)
    error('Some reference channels were not found in the EEG data.');
end

%% Compute Averaged Reference Channel
num_seconds = length(signal_data); % Total duration in seconds
reference_signal = cell(num_seconds, 1); % Initialize as a cell array

for j = 1:num_seconds
    temp_sum = zeros(size(signal_data{j})); % Initialize sum array for this second

    for idx = reference_indices
        temp_sum = temp_sum + signal_data{j}; % Sum the signals from reference channels
    end

    reference_signal{j} = temp_sum / length(reference_indices); % Compute average for this second
end

reference_signal = cell2mat(reference_signal);
aligned_signal_data = aligned_signal_data - reference_signal(1:length(aligned_signal_data));

%% Apply Band-Pass Filtering (0.1 - 40 Hz, Notch at 60 Hz)

low_cutoff = 0.1;  % Lower cutoff frequency
high_cutoff = 40;   % Upper cutoff frequency
order = 4;  % Filter order (higher order gives a sharper transition)
nyquist = sample_rate / 2;  
wn = [low_cutoff, high_cutoff] / nyquist;  % Normalized cutoff frequencies

% Design a band-pass Butterworth filter
[b_band, a_band] = butter(order, wn, 'bandpass');

% Apply zero-phase filtering to avoid phase distortion
bandpass_filtered_signal = filtfilt(b_band, a_band, aligned_signal_data);

% Manual Notch Filter Design (60 Hz)
notch_freq = 60;  % Frequency to be notched
bw = 2;  % Bandwidth of the notch filter

wo = notch_freq / nyquist;  % Normalized notch frequency
bw_norm = bw / nyquist;  % Normalized bandwidth

% Compute notch filter coefficients using second-order IIR formula
b_notch = [1, -2*cos(2*pi*wo), 1];
a_notch = [1, -2*(1 - bw_norm)*cos(2*pi*wo), (1 - bw_norm)^2];

% Apply the notch filter
filtered_signal_data = filtfilt(b_notch, a_notch, bandpass_filtered_signal);

%% Plot Frequency Response to Verify Filtering
[h_band, f_band] = freqz(b_band, a_band, 1024, sample_rate);
[h_notch, f_notch] = freqz(b_notch, a_notch, 1024, sample_rate);

figure;
subplot(2,1,1);
plot(f_band, abs(h_band));
xlabel('Frequency (Hz)');
ylabel('Gain');
title('Band-Pass Filter Frequency Response (0.1 - 40 Hz)');

subplot(2,1,2);
plot(f_notch, abs(h_notch));
xlabel('Frequency (Hz)');
ylabel('Gain');
title('Notch Filter Frequency Response (60 Hz)');

%% Optional: Plot Spectrogram to Check EEG Components
figure;
spectrogram(filtered_signal_data, 256, 250, 512, sample_rate, 'yaxis');
title('Filtered EEG Spectrogram (0.1 - 40 Hz and Notch Filter at 60 Hz)');

%% Ensure Time Vector Matches Signals (Prevents Plot Errors)

new_length = min([length(aligned_time_data), length(filtered_signal_data)]);
aligned_time_data = aligned_time_data(1:new_length);
aligned_signal_data = aligned_signal_data(1:new_length);
filtered_signal_data = filtered_signal_data(1:new_length);

%% Plot Filtered and Smoothed Signals

figure;
subplot(3,1,1);
plot(aligned_time_data, aligned_signal_data);
ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
title(['Original EEG Signal - ', full_signal_labels{channel_idx}]);
ylim([-1000 1000])
grid on;

subplot(3,1,2);
plot(aligned_time_data, filtered_signal_data, 'r');
ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
title(['Filtered EEG Signal - ', full_signal_labels{channel_idx}]);
ylim([-400 400])
grid on;


