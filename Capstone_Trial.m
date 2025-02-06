clc; close all; clear all;

% Authors: Kat

edf_file = 'C:\Users\sfmdf\OneDrive\Documents\MATLAB\Capstone\EEG1.edf'; % EDF file path
[hdr, record] = edfread(edf_file);
duration_minutes = input('Enter the duration in minutes to plot (1-60): ');

% Extract the signal labels (assuming EEG signals start from the 3rd column)
signal_labels = hdr.Properties.VariableNames(3:end);  % Adjust as needed

% Extract the EEG signal data (assuming the data for the first two channels are in columns 3 and 4)
signal_data = [hdr{:, 2}, hdr{:, 3}];  % Extract signals for C3 and C4

% Estimate sample rate (assumed to be 512 Hz)
sample_rate = 512; 
duration = height(record) / sample_rate;  % Total duration of the data in seconds

% Define time limit (first 2 minutes)
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
    
    % Create time for the current second (duplicate each second's time 512 times)
    second_time = time_per_second + (j-1);  % Align time with second j
    
    % Append to the aligned signal data and time
    aligned_signal_data = [aligned_signal_data; second_data_c3, second_data_c4];
    aligned_time_data = [aligned_time_data; second_time'];
end

% Plot the aligned signal data for both channels
figure;
for i = 1:2
    subplot(2, 1, i);
    plot(aligned_time_data, aligned_signal_data(:, i));
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['EEG Signal - ', signal_labels{i}]);
    grid on;
end
xlabel('Time (s)');