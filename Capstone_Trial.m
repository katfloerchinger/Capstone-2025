clc; close all; clear all;

% Authors: Kat

edf_file = 'C:\Users\sfmdf\OneDrive\Documents\MATLAB\Capstone\EEG1.edf'; % EDF file path
[hdr, record] = edfread(edf_file);
duration_minutes = height(hdr)/60;

% Extract the signal labels (assuming EEG signals start from the 3rd column)
signal_labels = hdr.Properties.VariableNames(3:end);  % Adjust as needed

% Extract the EEG signal data (assuming the data for the first two channels are in columns 3 and 4)
signal_data = [hdr{:, 2}, hdr{:, 3}];  % Extract signals for C3 and C4

% Estimate sample rate (assumed to be 512 Hz)
sample_rate = 512; 
time_limit = duration_minutes * 60; 
sample_limit = round(time_limit * sample_rate);  % Sample limit for 2 minutes

% Create a time vector for plotting each second (in terms of 512 points per second)
time_per_second = (0:511) / sample_rate;  % Time vector for 512 points in each second

% Create a new signal_data aligned with time for plotting
aligned_signal_data = [];
aligned_time_data = [];


% Loop through the data to align time with each 512 samples for every second
for j = 1:time_limit
    % Extract the 512 data points for the current second from the cell array
    second_data_c3 = signal_data{j, 1};  % Data for C3
    second_data_c4 = signal_data{j, 2};  % Data for C4
    second_time = time_per_second + (j-1);  % Align time with second j
    % Append to the aligned signal data and time
    aligned_signal_data = [aligned_signal_data; second_data_c3, second_data_c4];
    aligned_time_data = [aligned_time_data; second_time'];
end

% Plot the aligned signal data for both channels
figure;
for i = 1:2
    subplot(2, 1, i);
    plot(aligned_time_data(92100:1105920), aligned_signal_data(92100:1105920, i)); %plots from 3 minutes in (92100 samples aka 3 * 60 *512) until the end 
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['EEG Signal - ', signal_labels{i}]);
    grid on;
end
xlabel('Time (s)');


cutoff_freq = 0.5;  % Cutoff frequency in Hz
order = 4; % Filter order
[b, a] = butter(order, cutoff_freq/(sample_rate/2), 'low');
filtered_signal_data = filtfilt(b, a, aligned_signal_data);

for i = 1:2
    subplot(2, 2, i);
    plot(aligned_time_data(92100:1105920), aligned_signal_data(92100:1105920, i));
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Original EEG Signal - ', signal_labels{i}]);
    xlim([180 36*60])
    grid on;

    subplot(2, 2, i + 2);
    plot(aligned_time_data(92100:1105920), filtered_signal_data(92100:1105920, i), 'r');
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Filtered EEG Signal - ', signal_labels{i}]);
    xlim([180 36*60])
    grid on;
end
xlabel('Time (s)');


% Define moving average window size (in samples)
window_size = 512 * 5;  % 5-second window (512 samples per second)

% Compute running average (moving average)
smooth_signal_data(:, 1) = movmean(filtered_signal_data(92100:1105920, 1), window_size);
smooth_signal_data(:, 2) = movmean(filtered_signal_data(92100:1105920, 2), window_size);

for i = 1:2
    subplot(3, 2, i);
    plot(aligned_time_data(92100:1105920), aligned_signal_data(92100:1105920, i));
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Original EEG Signal - ', signal_labels{i}]);
    xlim([180 36*60])
    grid on;

    subplot(3, 2, i + 2);
    plot(aligned_time_data(92100:1105920), filtered_signal_data(92100:1105920, i), 'r');
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Filtered EEG Signal - ', signal_labels{i}]);
    xlim([180 36*60])
    grid on;

    subplot(3, 2, i + 4);
    plot(aligned_time_data(92100:1105920), smooth_signal_data(:, i), 'r');
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Smoothed EEG Signal - ', signal_labels{i}]);
    xlim([180 36*60])
    grid on;
end
xlabel('Time (s)');

% Phase space reconstruction parameters
tau = 1; % Time delay (arbitrary choice, can be optimized)
for i = 1:2
    figure;  % Create a new figure for each channel
    % Create the phase space plot
    plot(aligned_signal_data(2:end, i), aligned_signal_data(1:end-1, i), 'r');
    xlabel('Smooth Signal (t)', 'Interpreter', 'none');
    ylabel('Smooth Signal (t-1)', 'Interpreter', 'none');
    title(['Phase Space Plot for Channel ', signal_labels{i}], 'Interpreter', 'none');
    grid on;
end
title('Phase Space Reconstruction of EEG Data');

tau = 50;
figure;
aligned_signal_data = aligned_signal_data';
for ch = 1:2
    subplot(1,2,ch);
    X = aligned_signal_data(ch, :);
    X1 = X(1:end - 2*tau);
    X2 = X(1+tau:end - tau);
    X3 = X(1+2*tau:end);
    
    plot3(X1, X2, X3, 'b');
    grid on;
    xlabel('X(t)'); ylabel(['X(t+', num2str(tau), ')']); zlabel(['X(t+', num2str(2*tau), ')']);
    title(['EEG Channel ', num2str(ch)]);
end

sgtitle('Phase Space Reconstruction of EEG Data');