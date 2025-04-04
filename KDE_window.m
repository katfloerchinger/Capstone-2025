clc; close all; clear all;

%% Load EEG File

edf_file = '/Users/anthonytellez/Desktop/BIOEN 404/eeglab2024.2/EEG2.edf'; % EEG file path
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

%% Prompt User for Omission Time Bounds (Data Truncation)

start_omit_time = input('Enter how many seconds of data to omit from the beginning: ');
end_omit_time = input('Enter the end time (in seconds) after which data should be omitted: ');

% Validate inputs
if start_omit_time < 0 || end_omit_time <= start_omit_time || end_omit_time > aligned_time_data(end)
    error('Invalid time range. Ensure 0 <= start < end <= total duration.');
end

% Convert to sample indices
start_sample = round(start_omit_time * sample_rate) + 1;
end_sample = round(end_omit_time * sample_rate);

% Truncate all relevant arrays
aligned_signal_data = aligned_signal_data(start_sample:end_sample);
aligned_time_data = aligned_time_data(start_sample:end_sample);
aligned_ref_data = aligned_ref_data(start_sample:end_sample, :);

%% Re-referencing (on trimmed data)

ref = mean(aligned_ref_data(:, 1:end), 2);
re_ref_signal = aligned_signal_data - ref;

%% Apply Band-Pass Filtering (0.5 - 40 Hz)

low_cutoff = 0.5;  
high_cutoff = 40;  
order = 4;
nyquist = sample_rate / 2;
wn = [low_cutoff, high_cutoff] / nyquist;

[b, a] = butter(order, wn, 'bandpass');
filtered_signal_data = filtfilt(b, a, aligned_signal_data);

%% Compute Running Average (Moving Average)

window_size = 512 * 5;
filtered_signal = movmean(filtered_signal_data, window_size);

%% Phase Space Reconstruction

tau = 512;

X1 = filtered_signal(1+tau:end - tau);
X2 = filtered_signal(1+2*tau:end);
X3 = filtered_signal(1:end - 2*tau);

phase_data = [X1, X2, X3];
min_len = length(X1);
recon_time = aligned_time_data(1:min_len);  % Align time vector


%% KDE Analysis with Multiple Window Sizes and Weighted Averaging

% Define all window sizes (in seconds)
window_sizes_sec = [0.5, 1, 2, 5, 10];
master_window_sec = 20;

% Convert to samples
window_sizes_samples = window_sizes_sec * sample_rate;
master_window_samples = master_window_sec * sample_rate;

% Phase space data and time (from earlier)
usable_len = size(phase_data, 1);
recon_time = recon_time(1:usable_len);

% Storage for all KDE results
kde_results = struct();

for w = 1:length(window_sizes_sec)
    ws = window_sizes_sec(w);
    ws_samples = window_sizes_samples(w);
    num_windows = floor(usable_len / ws_samples);

    kde_vals = zeros(num_windows, 1);
    kde_times = zeros(num_windows, 1);

    for i = 1:num_windows
        idx_start = (i-1)*ws_samples + 1;
        idx_end = i*ws_samples;

        window_data = phase_data(idx_start:idx_end, :);
        norms = vecnorm(window_data, 2, 2);

        [f, ~] = ksdensity(norms);
        kde_vals(i) = max(f);
        kde_times(i) = recon_time(idx_end);
    end

    % Store results
    kde_results(w).kde_vals = kde_vals;
    kde_results(w).kde_times = kde_times;
    kde_results(w).window_size = ws;
end

%% Compute Final Aggregated KDE Values Using 20s Master Windows

% Determine master window bounds
num_master_windows = floor(usable_len / master_window_samples);
master_times = zeros(num_master_windows, 1);
final_kde_values = zeros(num_master_windows, 1);

for i = 1:num_master_windows
    m_start = (i-1)*master_window_samples + 1;
    m_end = i*master_window_samples;
    master_time = recon_time(m_end);
    master_times(i) = master_time;

    medians = [];

    for w = 1:length(window_sizes_sec)
        kde_vals = kde_results(w).kde_vals;
        kde_times = kde_results(w).kde_times;

        % Find values in this master window
        in_window = kde_times >= recon_time(m_start) & kde_times <= recon_time(m_end);

        if any(in_window)
            medians(end+1) = median(kde_vals(in_window));
        end
    end

    % Compute average of medians
    if ~isempty(medians)
        final_kde_values(i) = mean(medians);
    else
        final_kde_values(i) = NaN;
    end
end

%% Plot Final Weighted KDE Values Over Time (20s Master Windows)

figure;
plot(master_times, final_kde_values, '-o', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Weighted KDE Value');
title('Weighted KDE over Time (20s Windows from Multiple Resolutions)');
grid on;
