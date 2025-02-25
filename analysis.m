%% Enhanced Phase Space Analysis with Maximum Time Gap for Seizure Detection

% Setting parameters
filtered_Fz = filtered_Fz'; % reformating data
tolerance = 0.02;  % Tolerance for severity markers
test_window = 2000;  % Number of samples per analysis window
step_specificity = 1000;  % Step size between windows
tau = 512;  % Time delay for phase space reconstruction
alpha = 0.4; % compactness score weight
beta = 0.6; % linearization score weight

% Seizure Severity Marks
marks = [0.35, 0.45, 0.55, 0.7, 0.8];
messages = {
    'Predicting seizure activity.', ...
    'Possible progression into seizure severity 1.', ...
    'Possible progression into seizure severity 2.', ...
    'Probable seizure.', ...
    'Definitive seizure activity detected. Immediate medical attention required.'
};

% New: Define maximum allowed time gap for level 5 detections
max_gap_seconds = 5;  % If another level 5 marker appears within x secs, trigger emergency
last_marker5_time = -inf;  % Initialize last level 5 detection time to a very negative number

% Initialize tracking lists
compactness_score = []; 
linearization_score = []; 
combined_score = []; 
time_stamps = []; % Store corresponding time indices

% Score index for storing values properly
score_idx = 1;  

% Loop through EEG data in steps
for i = 1:step_specificity:(length(aligned_signal_data) - test_window - 2*tau)
    
    % make sure analysis is within index of data, break and notify
    if (i + test_window + 2*tau) > length(filtered_Fz)
        fprintf('Analysis for seizure goes over current known signal data.');
        break;  % Exit the loop to prevent out-of-bounds access
    end
    
    % Extract 3D Phase Space segment
    X1 = filtered_Fz(i:i + test_window, 1);  
    X2 = filtered_Fz(i+tau:i + tau + test_window, 1);  
    X3 = filtered_Fz(i+2*tau:i + 2*tau + test_window, 1);
    
    % Skip poorly conditioned segments
    if var(X1) < 1e-6 || var(X2) < 1e-6 || var(X3) < 1e-6
        continue;
    end

    % Compactness Test (Euclidean Distance)
    distances = sqrt((X1 - mean(X1)).^2 + (X2 - mean(X2)).^2 + (X3 - mean(X3)).^2);
    compactness = mean(distances);  % Lower distance → more compact phase space
    compactness = 1 - exp(-compactness);  % Normalize

    % Linearization Test (Least Squares Regression)
    % Center data for numerical stability
    X_centered = X1 - mean(X1);
    Y_centered = X2 - mean(X2);
    
    % Fit a linear model: Y = mX + b
    p = polyfit(X_centered, Y_centered, 1);  
    Y_fit = polyval(p, X_centered) + mean(X2);

    % Compute mean squared error (MSE)
    error_mse = mean((Y_centered - Y_fit).^2);
    linearization = 1 - exp(-error_mse);  % Normalize to 0-1 range

    % Combine Compactness + Linearization Scores
    combined = (compactness * alpha) + (linearization * beta);

    % Store values properly in the list
    compactness_score(score_idx) = compactness;
    linearization_score(score_idx) = linearization;
    combined_score(score_idx) = combined;
    time_stamps(score_idx) = i / sample_rate;  % Convert to seconds
    
    % Check Severity Thresholds
    for j = 1:length(marks)
        if combined >= marks(j) - tolerance
            fprintf('[ALERT] %s Time: %.2f seconds\n', messages{j}, time_stamps(score_idx));
        end
    end

    % New: Check Maximum Time Gap Between Level 5 Detections
    if combined >= marks(end) - tolerance  % If level 5 is detected
        current_time = time_stamps(score_idx);

        if (current_time - last_marker5_time) <= max_gap_seconds
            % If another level 5 detection is found within 10 sec, trigger emergency
            fprintf('[EMERGENCY] Patient experiencing ongoing seizure. Immediate medical attention needed.\n');
            fprintf('Time: %.2f seconds\n', current_time);
            break;
        end

        % Update last detected level 5 time
        last_marker5_time = current_time;
    end

    % Move to the next index
    score_idx = score_idx + 1;

end

% Plot Results
figure;
subplot(3,1,1);
plot(time_stamps, compactness_score, 'b');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Compactness Score');
title('Compactness Score Over Time');
grid on;

subplot(3,1,2);
plot(time_stamps, linearization_score, 'g');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Linearization Score');
title('Linearization Score Over Time');
grid on;

subplot(3,1,3);
plot(time_stamps, combined_score, 'm');
hold on;
yline(marks(end), 'r--', 'Seizure Threshold');
xlabel('Time (s)');
ylabel('Combined Score');
title('Combined Seizure Prediction Score');
grid on;
