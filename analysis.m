%% Analyzing phase space reconstruction

% setting which data to analyze
anyl_data = filtered_signal_data;

% Setting parameters
tolerance = 0.01;  % Linearization value tolerance 
test_window = 2000;  % Linearization analysis window
test_point_specificity = 40;  % Point gap
step_specificity = 1000;  % Gap step analysis

% Severity Marks
marks = [0.4, 0.5, 0.6, 0.75, 0.85];  
messages = {
    'Predicting seizure activity.', ...
    'Possible progression into seizure severity 1.', ...
    'Possible progression into seizure severity 2.', ...
    'Probable seizure.', ...
    'Definitive seizure activity detected. Immediate medical attention required.'
};

% Threshold for consecutive detection of definitive seizure activity
seizure_threshold = 20;  
seizure_counter = 0;  % Tracks consecutive detections

% Initialize linearization score tracking
linearization_score = zeros(length(anyl_data), 1);

% Loop through the data in chunks
for i = 1:step_specificity:length(anyl_data) - test_window
    % Extract segment of phase space
    x_segment = anyl_data(i:i + test_window, 1);  
    y_segment = anyl_data(i:i + test_window, 2);  

    % Skip poorly conditioned segments
    if var(x_segment) < 1e-6 || var(y_segment) < 1e-6
        continue;
    end

    % Center data to improve numerical stability
    x_centered = x_segment - mean(x_segment);
    y_centered = y_segment - mean(y_segment);

    % Fit linear model (robust to poorly conditioned data)
    p = polyfit(x_centered, y_centered, 1);  
    y_fit = polyval(p, x_centered) + mean(y_segment);
    
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
