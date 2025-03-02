% Video settings (MP4 format)
video_filename = 'phase_space_moving_window.mp4';
v = VideoWriter(video_filename, 'MPEG-4'); % Use MP4 format
v.FrameRate = 60; % Adjust playback speed
open(v);

% Prepare figure
figure;
hold on;
grid on;
xlabel('X(t)'); 
ylabel('X(t+\tau)'); 
zlabel('X(t+2\tau)');
title('Phase Space Evolution with Moving Window');

% Define parameters
tau = 512;  
window_seconds = 5; % Moving window duration
window_size = window_seconds * sample_rate; % Number of points in window

% Use entire filtered_Fz instead of a small segment
vid_segment = filtered_Fz(1:2300);
X = filtered_Fz';

X3 = X(1:end - 2*tau);
X2 = X(1+tau:end - tau);
X1 = X(1+2*tau:end);

num_points = length(X1);
h = plot3(NaN, NaN, NaN, 'b'); % Initialize empty plot

% Create video with moving window effect
for i = 1:num_points
    % Determine start of the rolling 20-second window
    start_idx = max(1, i - window_size);
    
    % Update plot with only the latest 20 seconds of data
    set(h, 'XData', X1(start_idx:i), ...
           'YData', X2(start_idx:i), ...
           'ZData', X3(start_idx:i));  
       
    drawnow;  % Update figure
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Finalize video
close(v);
disp('MP4 video saved successfully.');
