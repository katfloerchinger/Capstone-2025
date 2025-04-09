% Step 1: Load the EDF file
cfg = [];
cfg.dataset = 'EEG2.edf';
data = ft_preprocessing(cfg);

% OPTIONAL: Pick a channel to preview (e.g., channel 1)
channel_idx = 2;
time_full = cell2mat(data.time);
signal_full = data.trial{1}(channel_idx, :);

% Plot preview
figure;
plot(time_full, signal_full);
xlabel('Time (s)');
ylabel('Amplitude');
title(['Preview of EEG channel ', data.label{channel_idx}]);
grid on;

omit_time = input('Enter how many seconds of data to omit from the beginning: ');
omit_samples = omit_time * sample_rate;

if omit_samples < 0 || omit_samples >= length(aligned_signal_data)
    error('Invalid omission time. Must be within range of available data.');
end

omit_time_end = input('Enter how many seconds of data to omit from the end: ');
omit_samples_end = omit_time_end * sample_rate;

if omit_samples_end < 0 || omit_samples_end >= length(aligned_signal_data)
    error('Invalid omission time. Must be within range of available data.');
end

% Step 2: Define the trim window
start_time = omit_time_end;  % in seconds
end_time = data.time{1}(end) - omit_samples_end;

cfg = [];
cfg.latency = [start_time, end_time];
data_trimmed = ft_selectdata(cfg, data);

% Step 3: Adjust the time vector to preserve offset
time_offset = start_time;
data_trimmed.time{1} = data_trimmed.time{1} + time_offset;

% Optional: Update sampleinfo
num_samples = length(data_trimmed.time{1});
data_trimmed.sampleinfo = [time_offset * data.fsample + 1, ...
                           time_offset * data.fsample + num_samples];

% Step 4: Save trimmed data (EDF writing needs proper header)
cfg = [];
cfg.datafile = 'EEG2_Trim.edf';
cfg.headerfile = 'EEG2_Trim.vhdr';  % For other formats like BrainVision
cfg.dataformat = 'edf';