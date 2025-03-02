clc; close all; clear all;

edf_file = 'C:\Users\sfmdf\OneDrive\Documents\MATLAB\Capstone\EEG1.edf'; % EDF file path
[hdr, record] = edfread(edf_file);
duration_minutes = height(hdr) / 60;

% Extract the signal labels (assuming EEG signals start from the 3rd column)
signal_labels = hdr.Properties.VariableNames(3:end);  

% Extract EEG signal data

signal_data = [hdr{:, 13}, hdr{:, 14}, hdr{:, 15}, hdr{:, 11}, hdr{:, 12}, hdr{:, 20}, hdr{:, 21}, hdr{:, 22}, hdr{:, 23}, hdr{:, 4}, hdr{:, 5}]; 

% Estimate sample rate
sample_rate = 512; 
time_limit = duration_minutes * 60; 
sample_limit = round(time_limit * sample_rate);

% Range of seconds to test
seconds_to_test = 10;
start = 300000;
finish = start + seconds_to_test*512;

% Time vector for plotting
time_per_second = (0:511) / sample_rate;  

% Align signal data with time
aligned_signal_data = [];
aligned_time_data = [];

for j = 1:time_limit
    % Extract 512 samples for the current second
    Fz_second = signal_data{j, 1}; 
    other_channels = cell2mat(signal_data(j, 2:end)); % Convert other channels to matrix
    
    second_time = time_per_second + (j - 1);  
    aligned_signal_data = [aligned_signal_data; Fz_second, other_channels];
    aligned_time_data = [aligned_time_data; second_time'];
end

% Re-referencing: Remove common average from Fz
ref = mean(aligned_signal_data(:, 2:end), 2); 