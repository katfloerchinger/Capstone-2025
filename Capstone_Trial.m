% Authors: Kat

% MATLAB code to display EEG data from an EDF file

function display_eeg(edf_file)
    % Read EEG data from EDF file
    [hdr, record] = edfread(edf_file);
    
    % Get the number of channels and labels
    num_channels = size(record, 1);
    signal_labels = hdr.label;
    sample_rate = hdr.samples(1); % Assuming constant sample rate for all channels
    duration = size(record, 2) / sample_rate;
    
    % Time vector
    time = linspace(0, duration, size(record, 2));
    
    % Plot EEG signals
    figure;
    for i = 1:num_channels
        subplot(num_channels, 1, i);
        plot(time, record(i, :));
        ylabel(signal_labels{i}, 'Interpreter', 'none');
        title(['EEG Signal - ', signal_labels{i}]);
        grid on;
        if i < num_channels
            set(gca, 'XTickLabel', []);
        end
    end
    xlabel('Time (s)');
    suptitle('EEG Signals');
end

% Example usage
edf_file = 'C:\Users\sfmdf\OneDrive\Documents\MATLAB\EEG1.edf'; % EDF file path
display_eeg(edf_file);
