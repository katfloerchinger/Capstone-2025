%% **Low-Pass Filter Design**
cutoff_freq = 0.5;  % Cutoff frequency in Hz
order = 4; % Filter order

% Design Butterworth low-pass filter
[b, a] = butter(order, cutoff_freq/(sample_rate/2), 'low');

% Apply the filter to each channel
filtered_signal_data = filtfilt(b, a, aligned_signal_data);

%% **Plot Original and Filtered EEG Signals**
figure;

for i = 1:2
    subplot(2, 2, i);
    plot(aligned_time_data, aligned_signal_data(:, i));
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Original EEG Signal - ', signal_labels{i}]);
    grid on;

    subplot(2, 2, i + 2);
    plot(aligned_time_data, filtered_signal_data(:, i), 'r');
    ylabel(signal_labels{i}, 'Interpreter', 'none');
    title(['Filtered EEG Signal - ', signal_labels{i}]);
    grid on;
end

xlabel('Time (s)');
