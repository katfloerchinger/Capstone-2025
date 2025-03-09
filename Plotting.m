close all

% data_type_for_plot = smooth_signal_data;
%data_type_for_plot = aligned_signal_data;
data_type_for_plot = filtered_signal_data;

%% Phase Space Reconstruction in 45-Second Increments (1100s - 1280s)

tau = 512; % Time delay (adjust if needed)
increment = 5;
increment_duration = increment * sample_rate; % 45-second segments
start_t = 1250;
start_time = start_t * sample_rate; % Convert seconds to samples
end_t = 1300;
end_time = end_t * sample_rate; % Convert seconds to samples
num_segments = floor((end_time - start_time) / increment_duration);

figure(1);
for i = 1:num_segments
    % Extract segment
    segment_start = start_time + (i - 1) * increment_duration;
    segment_end = segment_start + increment_duration - 1;
    segment_data = data_type_for_plot(segment_start:segment_end);
    
    % Create phase space vectors
    X = segment_data';
    if length(X) < 2*tau
        continue; % Skip if too short
    end
    X3 = X(1:end - 2*tau); % Backward
    X2 = X(1+2*tau:end);   % Forward
    X1 = X(1+tau:end - tau); % Standard

    % Subplot for each 45-second segment
    subplot(ceil(num_segments / 2), 2, i); % Organizes in 2 columns
    plot3(X1, X2, X3, 'b');
    grid on;
    xlabel('X(t)');
    ylabel('X(t + 1)'); %ylabel(['X(t+', num2str(tau), ')']);
    zlabel('X(t + 2)'); %zlabel(['X(t+', num2str(2*tau), ')']);
    title(['Phase Space (', num2str(start_t + (i - 1) * increment), '-', num2str(end_t + i * increment), 's)']);
end
sgtitle(['Phase Space Reconstruction (' num2str(start_t) ' - ' num2str(end_t) ' s) - ' full_signal_labels{channel_idx}]);


%% Define the Time Range for 10-Second Segment
start_time = 1241; 
end_time = 1245;

start_idx = start_time * sample_rate + 1;
end_idx = end_time * sample_rate;
time_segment = aligned_time_data(start_idx:end_idx);
signal_segment = data_type_for_plot(start_idx:end_idx);

figure;
subplot(4,1,1); % Takes 1/4th of the figure
plot(time_segment, signal_segment, 'b');
xlabel('Time (s)');
ylabel('EEG');
title(['EEG Signal (', num2str(start_time), 's - ', num2str(end_time), 's)']);
grid on;
set(gca, 'XTick', [], 'XColor', 'none'); % Hide X-axis ticks

subplot(4,1,[2 3 4]); % Takes 3/4th of the figure
tau = 512;  

X = signal_segment'; 
X3 = X(1:end - 2*tau);
X2 = X(1+tau:end - tau);
X1 = X(1+2*tau:end);

plot3(X1, X2, X3, 'r');
grid on;
xlabel('X(t)');
ylabel('X(t+\tau)');
zlabel('X(t+2\tau)');
title(['Phase Space Reconstruction- ' full_signal_labels{channel_idx}]);
