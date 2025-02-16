clc; clear; close all;

% Generate sample EEG data (5 channels, 1000 samples each)
num_channels = 5;
num_samples = 1000;
t = linspace(0, 10, num_samples); % Simulated time vector (10 seconds)
fs = num_samples / 10; % Sampling frequency (100 Hz)

% Create synthetic EEG data (sine waves with noise)
eeg_data = sin(2 * pi * (1:num_channels)' * t) + 0.1 * randn(num_channels, num_samples);
tau = 10; 


% Plot phase space reconstruction for each channel
figure;
for ch = 1:num_channels
    subplot(2,3,ch);
    X = eeg_data(ch, :);
    X1 = X(1:end - 2*tau);
    X2 = X(1+tau:end - tau);
    X3 = X(1+2*tau:end);
    
    plot3(X1, X2, X3, 'b');
    grid on;
    xlabel('X(t)'); ylabel(['X(t+', num2str(tau), ')']); zlabel(['X(t+', num2str(2*tau), ')']);
    title(['EEG Channel ', num2str(ch)]);
end

sgtitle('Phase Space Reconstruction of EEG Data');
