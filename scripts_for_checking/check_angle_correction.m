% Define the sequence of orientations including jumps
testOrientations = [10, 20, 30, 40, 50, 80, -85, -75, -65, -55, -45, -10, ...
                    1, 10, 30, 50, 70, 80, 87, -87, -86,-60,-50,-30,-20,-10, ...
                    10,30,10,3,-3,-7,-20,-60,-80,85,80,60,30,20,-10,-20,-20,-30,-40,-50,-50,-50,-40,-30,-20,-10,10,10,10,10,10,10,10];
numFrames = length(testOrientations);

% Initialize the array for corrected cumulative angles
cumulativeAngles = zeros(1, numFrames);

% Initialize cumulative angle with the first orientation
cumulativeAngles(1) = testOrientations(1);

% Process each frame to correct orientations
for frame = 2:numFrames
    % Calculate the angle difference
    angleDifference = testOrientations(frame) - testOrientations(frame - 1);
    
    % Correct for flipping
    if abs(angleDifference) > 90
        if angleDifference > 0
            angleDifference = angleDifference - 180;
        else
            angleDifference = angleDifference + 180;
        end
    end
    
    % Update the cumulative angle
    cumulativeAngles(frame) = cumulativeAngles(frame - 1) + angleDifference;
end

% Plot the original and cumulative corrected orientations
figure;
subplot(2,1,1);
plot(1:numFrames, testOrientations, 'r-', 'LineWidth', 2, 'DisplayName', 'Original Orientations');
xlabel('Frame Number');
ylabel('Orientation (degrees)');
legend('show');
title('Original Orientations');
grid on;

subplot(2,1,2);
plot(1:numFrames, cumulativeAngles, 'b-', 'LineWidth', 2, 'DisplayName', 'Cumulative Corrected Angles');
xlabel('Frame Number');
ylabel('Cumulative Angle (degrees)');
legend('show');
title('Cumulative Corrected Angles');
grid on;
