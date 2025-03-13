clc; close all; clear all;

%% Load EEG File

edf_file = '/Users/anthonytellez/Desktop/BIOEN 404/eeglab2024.2/EEG2.edf'; % EEG file path
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

%% Ask User for Data Retention Time Range
retain_start_time = input('Enter the start time (in seconds) to KEEP the data from: ');
retain_end_time = input('Enter the end time (in seconds) to KEEP the data until: ');

% Convert to number of data points
retain_start_samples = round(retain_start_time * sample_rate);
retain_end_samples = round(retain_end_time * sample_rate);

% Ensure valid input
if retain_start_samples < 1 || retain_end_samples > length(aligned_signal_data) || retain_start_samples >= retain_end_samples
    error('Invalid time range. Ensure start is before end and within the available data.');
end

% Keep only the selected range of data
aligned_signal_data = aligned_signal_data(retain_start_samples:retain_end_samples);
aligned_time_data = aligned_time_data(retain_start_samples:retain_end_samples);

% for EEG 2, start at 520, end at 2380

%% Apply Band-Pass Filtering (0.1 - 40 Hz)

low_cutoff = 0.1;  % per Dr Shahin
high_cutoff = 40;  % need to mess around with
order = 4;  % Filter order
nyquist = sample_rate / 2;  
wn = [low_cutoff, high_cutoff] / nyquist;  % Normalized cutoff frequencies

% Design a band-pass Butterworth filter
[b, a] = butter(order, wn, 'bandpass');

% Apply zero-phase filtering
filtered_signal_data = filtfilt(b, a, aligned_signal_data);

%% Ensure Time Vector Matches Signals (Prevents Plot Errors)

new_length = min([length(aligned_time_data), length(filtered_signal_data)]);
aligned_time_data = aligned_time_data(1:new_length);
aligned_signal_data = aligned_signal_data(1:new_length);
filtered_signal_data = filtered_signal_data(1:new_length);

%% Multi-Window Phase Space Reconstruction with Multiple Point Density Methods

while true
    % Ask user for tau
    tau_input = input('Enter tau (time delay in samples, e.g., 512) or type "end" to stop: ', 's');
    
    % Exit condition
    if strcmpi(tau_input, "end")
        disp('Exiting Phase Space Visualization.');
        break;
    end
    
    % Convert tau input to a number
    tau = str2double(tau_input);
    
    % Ensure tau is valid
    if isnan(tau) || tau <= 0 || tau >= length(filtered_signal_data)
        disp('Invalid tau value. Tau must be a positive number within the data range.');
        continue;
    end

    % Ask the user for the analysis window range
    range_seconds = input('Enter the time window range in seconds (e.g., 20): ');

    % Ensure valid input
    if isempty(range_seconds) || range_seconds <= 0
        disp('Invalid range. Please enter a positive number.');
        continue;
    end

    % Ask the user for the central time point
    center_time = input('Enter the central time in seconds (e.g., 1300): ') - retain_start_time;

    % Convert time to sample indices
    center_idx = round(center_time * sample_rate);
    range_samples = round(range_seconds * sample_rate);

    % Define 3 pre-center and 4 post-center time windows
    time_windows = [-3, -2, -1, 0, 1, 2, 3, 4] * range_samples + center_idx;

    % Initialize data storage
    all_phase_space_points = [];  
    density_bounding_box = [];
    density_convex_hull = [];
    density_kde = [];

    % Create figure for visualization
    figure;
    sgtitle('Multi-Window Phase Space Reconstruction');

    plot_index = 1;
    for i = 1:length(time_windows) - 1
        start_idx = time_windows(i);
        end_idx = start_idx + range_samples;
        
        % Ensure the time range is within bounds
        if start_idx < 1 || end_idx > length(filtered_signal_data)
            fprintf('Skipping window: %d to %d (out of bounds)\n', start_idx, end_idx);
            continue; 
        end
        
        % Extract the selected time range
        PSR_data = filtered_signal_data(start_idx:end_idx)';

        % Compute delayed vectors
        X1 = PSR_data(1+tau:end - tau);
        X2 = PSR_data(1+2*tau:end);
        X3 = PSR_data(1:end - 2*tau);

        % Store phase space points in column format
        phase_space_points = [X1; X2; X3];

        % Append phase space points to the full dataset
        all_phase_space_points = [all_phase_space_points, phase_space_points];

        % Compute point densities
        density_bounding_box(end+1) = computeBoundingBoxDensity(phase_space_points);
        density_convex_hull(end+1) = computeConvexHullDensity(phase_space_points);
        density_kde(end+1) = computeKDEDensity(phase_space_points);

        % Debugging: Print densities
        fprintf('Window: %d to %d, Bounding Box: %f, Convex Hull: %f, KDE: %f\n', ...
            start_idx, end_idx, density_bounding_box(end), density_convex_hull(end), density_kde(end));

        % Plot phase space
        subplot(2, 4, plot_index);
        plot3(X1, X2, X3, 'r', 'LineWidth', 1);
        grid on;
        xlabel('X(t)');
        ylabel(['X(t+', num2str(tau), ')']);
        zlabel(['X(t-', num2str(tau), ')']);
        title(['Range: ', num2str((start_idx/sample_rate) + retain_start_time), 's - ', num2str((end_idx/sample_rate) + retain_start_time), 's']);

        plot_index = plot_index + 1;
    end

    % Plot all density methods over time
    time_labels = (time_windows(1:end-1) / sample_rate) + retain_start_time;

    figure;
    hold on;
    plot(time_labels, density_bounding_box, '-o', 'DisplayName', 'Bounding Box Density');
    plot(time_labels, density_convex_hull, '-o', 'DisplayName', 'Convex Hull Density');
    plot(time_labels, density_kde, '-o', 'DisplayName', 'KDE Density');
    xlabel('Time Window (s)');
    ylabel('Point Density');
    title('Comparison of Point Density Methods Over Time');
    legend;
    grid on;
end

%% Bounding Box Density Calculation
function density = computeBoundingBoxDensity(points)
    if size(points, 2) < 4
        density = NaN;
        return;
    end

    x_min = min(points(1,:)); x_max = max(points(1,:));
    y_min = min(points(2,:)); y_max = max(points(2,:));
    z_min = min(points(3,:)); z_max = max(points(3,:));

    volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    if volume <= 0
        density = NaN;
    else
        density = size(points, 2) / volume;
    end
end

%% Convex Hull Density Calculation
function density = computeConvexHullDensity(points)
    if size(points, 2) < 4
        density = NaN;
        return;
    end
    try
        [~, volume] = convhull(points');
        density = size(points, 2) / volume;
    catch
        density = NaN;  % Handle convex hull failures
    end
end

%% KDE Density Calculation
function density = computeKDEDensity(points)
    if size(points, 2) < 4
        density = NaN;
        return;
    end
    try
        [~, density_vals] = ksdensity(points(:));  % Flatten for KDE
        density = mean(density_vals);
    catch
        density = NaN;
    end
end
