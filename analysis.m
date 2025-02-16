%% Analyzing phase space reconstruction

% Setting parameters
tolerance = 0.01;  % Linearization value tolerance 
test_window = 100;  % Linearization analysis window
test_point_specificity = 1;  % Point gap
step_specificity = 10;  % Gap step analysis

% Severity Marks
marks = [0.4, 0.5, 0.6, 0.75, 0.85];  
messages = {
    'Predicting seizure activity.', ...
    'Possible progression into seizure severity 1.', ...
    'Possible progression into seizure severity 2.', ...
    'Probably seizure.', ...
    'Definitive seizure activity detected. Immediate medical attention required.'
};

% Threshold for consecutive detection of definitive seizure activity
seizure_threshold = 1000;  
seizure_counter = 0;  % Tracks consecutive detections

% Initialize linearization score tracking
linearization_score = zeros(length(aligned_signal_data), 1);

% Loop through the data in chunks
for i = 1:step_specificity:length(aligned_signal_data) - test_window
    % Extract the segment of phase space
    x_segment = aligned_signal_data(i:i + test_window, 1);  % Channel C3
    y_segment = aligned_signal_data(i:i + test_window, 2);  % Channel C4

    % Fit a linear model (y = ax + b)
    p = polyfit(x_segment, y_segment, 1);  
    y_fit = polyval(p, x_segment);
    
    % Compute the error (how linear the phase space is)
    error_mse = mean((y_segment - y_fit).^2);
    
    % Normalize the error (scaling to 0-1 range)
    linearization_score(i) = 1 - exp(-error_mse); 

    % Check if score exceeds severity thresholds
    for j = 1:length(marks)
        if linearization_score(i) >= marks(j) - tolerance
            fprintf('[ALERT] %s Time: %.2f seconds\n', messages{j}, i / sample_rate);
        end
    end

    % Emergency condition: Check for consecutive detections
    if linearization_score(i) >= marks(end) - tolerance
        seizure_counter = seizure_counter + 1;  % Increment if marker 5 is triggered
    else
        seizure_counter = 0;  % Reset if marker 5 is not triggered consecutively
    end

    % If marker 5 is detected consecutively, trigger emergency response and break
    if seizure_counter >= seizure_threshold
        fprintf('[EMERGENCY] Patient experiencing ongoing seizure. Immediate medical attention needed.\n');
        fprintf('Time: %.2f seconds\n', i / sample_rate);
        break;
    end
end
