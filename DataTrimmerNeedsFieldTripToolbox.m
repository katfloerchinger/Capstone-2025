% Step 1: Load the original EDF using FieldTrip
cfg = [];
cfg.dataset = 'EEG2.edf';
data = ft_preprocessing(cfg);

% Step 2: Manually trim the data
fs = data.fsample;
omit_samples_start = input('Seconds to omit at start: ') * fs;
omit_samples_end = input('Seconds to omit at end: ') * fs;

data_trimmed = data;  % Copy the structure
for i = 1:length(data.trial)
    data_trimmed.trial{i} = data.trial{i}(:, omit_samples_start+1:end - omit_samples_end);
    data_trimmed.time{i} = data.time{i}(omit_samples_start+1:end - omit_samples_end);
end

% Step 3: Reconstruct a basic header (FieldTrip format)
hdr = [];
hdr.Fs = fs;
hdr.nChans = size(data_trimmed.trial{1}, 1);
hdr.label = data_trimmed.label;
hdr.nSamples = size(data_trimmed.trial{1}, 2);
hdr.nTrials = 1;
hdr.label = data.label;
hdr.orig = [];  % Not required

% Step 4: Write trimmed data to EDF
ft_write_data('EEG2_Trim.edf', data_trimmed.trial{1}, 'header', hdr, 'dataformat', 'edf');
der', data_trimmed.hdr);
