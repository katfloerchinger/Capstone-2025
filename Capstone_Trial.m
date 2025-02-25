clc; close all; clear all;

% Authors: Kat/Joey

edf_file = '/Users/joeyroberts/Desktop/CAPSTONE/EEG2.edf'; % EDF file path
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
Fz = aligned_signal_data(:, 1) - ref;

% Plot raw Fz signal
figure;
plot(aligned_time_data(start:finish), Fz(start:finish)); 
ylabel('Fz', 'Interpreter', 'none');
title('EEG Signal - Fz');
grid on;
xlabel('Time (s)');

%% **High-Pass Filter (0.5 Hz)**
low_cutoff = 0.5;  
order = 4;  

[b, a] = butter(order, low_cutoff / (sample_rate / 2), 'high');
high_passed_Fz = filtfilt(b, a, Fz);

%% **Gaussian Low-Pass Filter (100 Hz)**
sigma = 10; 
window_size = 50; 

% Create a 1D Gaussian kernel
gaussian_kernel = gausswin(window_size);
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize

% Apply Gaussian filter using convolution
filtered_Fz = conv(high_passed_Fz, gaussian_kernel, 'same');

%% **Plot Filtered Fz Signal**
figure;
plot(aligned_time_data(start:finish), filtered_Fz(start:finish), 'r');
ylabel('Fz', 'Interpreter', 'none');
title('Filtered EEG Signal - Fz');
grid on;
xlabel('Time (s)');

%% **Phase Space Reconstruction for Fz**
tau = 256;
figure;
X = filtered_Fz(start:finish)';  
X3 = X(1:end - 2*tau);
X2 = X(1+tau:end - tau);
X1 = X(1+2*tau:end);

plot3(X1, X2, X3, 'b');
grid on;
xlabel('X(t)'); 
ylabel(['X(t+', num2str(tau), ')']); 
zlabel(['X(t+', num2str(2*tau), ')']);
title('Phase Space Reconstruction - Fz');