clc; close all; clear all;

%% Load EEG File

edf_file = '/Users/joeyroberts/Desktop/CAPSTONE/EEG2.edf'; % EEG file path
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

ref_data = [hdr{:, 14}, hdr{:, 15}, hdr{:, 11}, hdr{:, 12}, hdr{:, 20}, hdr{:, 21}, hdr{:, 22}, hdr{:, 23}, hdr{:, 4}, hdr{:, 5}]; 

%% Estimate Sample Rate and Create Time Vector

sample_rate = 512;
time_limit = duration_minutes * 60;
sample_limit = round(time_limit * sample_rate);
time_per_second = (0:511) / sample_rate;

%% Align Signal Data with Time

aligned_signal_data = [];
aligned_time_data = [];
aligned_ref_data = [];

for j = 1:time_limit
    second_data = signal_data{j};  
    ref_sec_data = cell2mat(ref_data(j, 1:end));
    second_time = time_per_second + (j-1);
    aligned_signal_data = [aligned_signal_data; second_data];
    aligned_ref_data = [aligned_ref_data; ref_sec_data];
    aligned_time_data = [aligned_time_data; second_time'];
end

%% Display Raw EEG Signal Before Omitting Data

figure;
plot(aligned_time_data, aligned_signal_data);
ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
title(['EEG Signal - ', full_signal_labels{channel_idx}]);
grid on;
xlabel('Time (s)');

omit_time = input('Enter how many seconds of data to omit from the beginning: ');
omit_samples = omit_time * sample_rate;

if omit_samples < 0 || omit_samples >= length(aligned_signal_data)
    error('Invalid omission time. Must be within range of available data.');
end

aligned_signal_data(1:omit_samples) = [];
aligned_time_data(1:omit_samples) = [];
aligned_ref_data(1:omit_samples) = [];

%% Re-referencing

ref = mean(aligned_ref_data(:, 1:end), 2);
re_ref_signal = aligned_signal_data - ref;

%% Apply Band-Pass Filtering (0.5 - 100 Hz)

low_cutoff = 0.5;  
high_cutoff = 100;  
order = 4;  % Filter order
nyquist = sample_rate / 2;  
wn = [low_cutoff, high_cutoff] / nyquist;  % Normalized cutoff frequencies

% Design a band-pass Butterworth filter
[b, a] = butter(order, wn, 'bandpass');

% Apply zero-phase filtering
filtered_signal_data = filtfilt(b, a, aligned_signal_data);

%% Compute Running Average (Moving Average)

window_size = 512 * 5;  % 5-second window
filtered_signal = movmean(filtered_signal_data, window_size);

%% Subplot Filtered Signal in 20-Second Intervals

seconds_to_test = 20;
num_intervals = 6;
start_times = linspace(650000, 650000 + (num_intervals - 1) * 20 * sample_rate, num_intervals);

figure;
for i = 1:num_intervals
    start = round(start_times(i));
    finish = start + seconds_to_test * sample_rate;
    subplot(2, 3, i);
    plot(aligned_time_data(start:finish), filtered_signal(start:finish), 'r');
    ylabel(full_signal_labels{channel_idx}, 'Interpreter', 'none');
    title(['Filtered EEG Signal - Interval ', num2str(i)]);
    grid on;
    xlabel('Time (s)');
end

%% Phase Space Reconstruction

tau = 256;
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
    grid on;
    xlabel('X(t)'); ylabel(['X(t-', num2str(tau/512), 'sec)']); zlabel(['X(t-', num2str(2*(tau/512)), 'sec)']);
    title(['Phase Space Reconstruction - Interval ', num2str(i)]);
end

%% Enhanced Phase Space Analysis with Maximum Time Gap for Seizure Detection

% Setting parameters
smooth_signal_data = filtered_signal'; % reformating data
tolerance = 0.02;  % Tolerance for severity markers
test_window = 2000;  % Number of samples per analysis window
step_specificity = 1000;  % Step size between windows
tau = 256;  % Time delay for phase space reconstruction
gamma = 2.5;  % Sensitivity to sudden compactness changes

% Seizure Severity Marks
marks = [0.35, 0.45, 0.55, 0.7, 0.8];
messages = {
    'Predicting seizure activity.', ...
    'Possible progression into seizure severity 1.', ...
    'Possible progression into seizure severity 2.', ...
    'Probable seizure.', ...
    'Definitive seizure activity detected. Immediate medical attention required.'
};

% New: Define maximum allowed time gap for level 5 detections
max_gap_seconds = 30;  % If another level 5 marker appears within x secs, trigger emergency
last_marker5_time = -inf;  % Initialize last level 5 detection time to a very negative number

% Initialize tracking lists
compactness_score = []; 
linearization_score = []; 
combined_score = []; 
time_stamps = []; % Store corresponding time indices

% Score index for storing values properly
score_idx = 1;  

% Loop through EEG data in steps
for i = 1:step_specificity:(length(aligned_signal_data) - test_window - 2*tau)
    
    % make sure analysis is within index of data, break and notify
    if (i + test_window + 2*tau) > length(smooth_signal_data)
        fprintf('Analysis for seizure goes over current known signal data.');
        break;  % Exit the loop to prevent out-of-bounds access
    end
    
    % Extract 3D Phase Space segment
    X1 = smooth_signal_data(i + tau:i + test_window - tau); % standard
    X2 = smooth_signal_data(i+ 2*tau:i + test_window); % forward
    X3 = smooth_signal_data(i:i - 2*tau + test_window); % back
    
    % Skip poorly conditioned segments
    if var(X1) < 1e-6 || var(X2) < 1e-6 || var(X3) < 1e-6
        continue;
    end

    % Compactness Test (Euclidean Distance)
    distances = sqrt((X1 - mean(X1)).^2 + (X2 - mean(X2)).^2 + (X3 - mean(X3)).^2);
    compactness = mean(distances);  % Lower distance → more compact phase space
    compactness = 1 - exp(-compactness);  % Normalize

    % Linearization Test (Least Squares Regression)
    % Center data for numerical stability
    X_centered = X1 - mean(X1);
    Y_centered = X2 - mean(X2);
    
    % Fit a linear model: Y = mX + b
    p = polyfit(X_centered, Y_centered, 1);  
    Y_fit = polyval(p, X_centered) + mean(X2);

    % Compute mean squared error (MSE)
    error_mse = mean((Y_centered - Y_fit).^2);
    linearization = 1 - exp(-error_mse);  % Normalize to 0-1 range

    % Compute weightings dynamically
    w_C = exp(-linearization);  % Compactness dominates when L is low
    w_L = 1 - exp(-compactness);  % Linearization dominates when C is high

    % Compute compactness change rate (handle first iteration)
    if score_idx > 1
        delta_C = abs(compactness - compactness_score(score_idx - 1));
    else
        delta_C = 0;
    end

    % Compute enhanced combined score with sigmoid normalization
    combined = 1 / (1 + exp(-(w_C * compactness + w_L * linearization + gamma * delta_C)));

    % Store values properly in the list
    compactness_score(score_idx) = compactness;
    linearization_score(score_idx) = linearization;
    combined_score(score_idx) = combined;
    time_stamps(score_idx) = i / sample_rate;  % Convert to seconds
    
    % Check Severity Thresholds
    for j = 1:length(marks)
        if combined >= marks(j) - tolerance
            fprintf('[ALERT] %s Time: %.2f seconds\n', messages{j}, time_stamps(score_idx)+omit_time);
        end
    end

    % New: Check Maximum Time Gap Between Level 5 Detections
    if combined >= marks(end) - tolerance  % If level 5 is detected
        current_time = time_stamps(score_idx);

        if (current_time - last_marker5_time) <= max_gap_seconds
            % If another level 5 detection is found within 10 sec, trigger emergency
            fprintf('[EMERGENCY] Patient experiencing ongoing seizure. Immediate medical attention needed.\n');
            fprintf('Time: %.2f seconds\n', current_time+omit_time);
            break;
        end

        % Update last detected level 5 time
        last_marker5_time = current_time;
    end

    % Move to the next index
    score_idx = score_idx + 1;

end

% Plot Results
figure;
subplot(3,1,1);
plot(time_stamps+omit_time, compactness_score, 'b');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Compactness Score');
title('Compactness Score Over Time');
grid on;

subplot(3,1,2);
plot(time_stamps+omit_time, linearization_score, 'g');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Linearization Score');
title('Linearization Score Over Time');
grid on;

subplot(3,1,3);
plot(time_stamps+omit_time, combined_score, 'm');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Combined Score');
title('Combined Seizure Prediction Score');
grid on;
