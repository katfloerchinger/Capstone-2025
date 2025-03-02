% Authors: Kat

% Re-referencing: Remove common average from Fz
Fz = aligned_signal_data(:, 1) - ref;

% Plot raw Fz signal
figure;
plot(aligned_time_data(start:finish), Fz(start:finish)); 
ylabel('Fz', 'Interpreter', 'none');
title('EEG Signal - Fz');
grid on;
xlabel('Time (s)');

%% **High-Pass Filter (0.4 Hz)**
low_cutoff = 0.4;  
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
tau = 512;
figure;
X = filtered_Fz(start:finish)';  
X3 = X(1:end - 2*tau);
X2 = X(1+tau:end - tau);
X1 = X(1+2*tau:end);

plot3(X1, X2, X3, 'b');
grid on;
xlabel('X(t)'); 
ylabel('X(t+1)'); 
zlabel('X(t+2)');
title('Phase Space Reconstruction - Fz');