% Define band-pass cutoff frequencies
low_cutoff = 0.5;   % Lower cutoff frequency (Hz)
high_cutoff = 100;  % Upper cutoff frequency (Hz)
order = 4;          % Filter order

% Normalize the cutoff frequencies (divide by Nyquist frequency)
nyquist = sample_rate / 2;  
wn = [low_cutoff, high_cutoff] / nyquist;  % Normalized cutoff frequencies

% Design a band-pass Butterworth filter
[b, a] = butter(order, wn, 'bandpass');

% Apply zero-phase filtering using filtfilt (prevents phase distortion)
filtered_signal_data = filtfilt(b, a, aligned_signal_data);
