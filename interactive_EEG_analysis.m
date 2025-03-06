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

%% Apply Band-Pass Filtering (0.5 - 1000 Hz)

low_cutoff = 0.1;  % per Dr Shahin
high_cutoff = 60;  % need to mess around with
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

%% Multi-Window EEG analysis
while true
    % Ask user to continue wuth another analysis
    cont_input = input('Input "end" to complete analysis. Hit "enter" to continue. ', 's');
    
    % Exit condition
    if strcmpi(cont_input, "end") | strcmpi(cont_input, "stop")
        disp('Exiting EEG Analysis.');
        break;
    end

    % Ask user for the analysis window range
    range_seconds = input('Enter the time window range in seconds (e.g., 20): ');

    % Ensure valid input
    if isempty(range_seconds) || range_seconds <= 0
        disp('Invalid range. Please enter a positive number.');
        continue;
    end

    % Ask user for the central time point
    center_time = input('Enter the central time in seconds (e.g., 1300): ');

    % Convert time to sample indices
    center_idx = round(center_time * sample_rate);
    range_samples = round(range_seconds * sample_rate);

    % Define 3 pre-center and 4 post-center time windows
    time_windows = [-3, -2, -1, 0, 1, 2, 3, 4] * range_samples + center_idx;

    % Initialize periodicity scores
    periodicity_scores = nan(1, length(time_windows) - 1);
    valid_time_stamps = nan(1, length(time_windows) - 1); % Store valid time windows

    % Create figure for visualization
    figure;
    sgtitle('Multi-Window EEG Periodicity Analysis');

    plot_index = 1; % Track subplot index
    valid_idx = 1; % Track valid periodicity scores

    for i = 1:length(time_windows)
        start_idx = time_windows(i);
        end_idx = start_idx + range_samples;

        % Ensure time range is within bounds
        if start_idx < 1 || end_idx > length(filtered_signal_data)
            continue; % Skip invalid windows
        end
        
        % Extract EEG segment for periodicity analysis
        EEG_segment = filtered_signal_data(start_idx:end_idx);

        %% Compute Periodicity Score (Mean Recurrence Rate)
        dist_matrix = squareform(pdist(EEG_segment', 'euclidean'));
        threshold = 0.1 * mean(dist_matrix(:));  
        recurrence_matrix = dist_matrix < threshold;
        periodicity_scores(valid_idx) = mean(recurrence_matrix(:));
        valid_time_stamps(valid_idx) = start_idx / sample_rate; % Store corresponding time

        % Plot EEG segment
        subplot(2, 4, plot_index);
        plot((0:length(EEG_segment)-1) / sample_rate, EEG_segment, 'r', 'LineWidth', 1);
        grid on;
        xlabel('Time (s)');
        ylabel('EEG Signal');
        title(['Range: ', num2str(start_idx/sample_rate), 's - ', num2str(end_idx/sample_rate), 's']);

        plot_index = plot_index + 1;
        valid_idx = valid_idx + 1;
    end

    % Remove NaN values from periodicity_scores and valid_time_stamps
    valid_indices = ~isnan(periodicity_scores);
    periodicity_scores = periodicity_scores(valid_indices);
    valid_time_stamps = valid_time_stamps(valid_indices);

    % Ensure we have at least one valid periodicity score before plotting
    if ~isempty(valid_time_stamps)
        figure;
        plot(valid_time_stamps, periodicity_scores, '-o', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Periodicity Score');
        title('Periodicity Score Across EEG Time Windows');
        grid on;
    else
        disp('No valid periodicity scores computed. Try a different time range.');
    end

    % Pause before next iteration
    pause(0.3);
end
