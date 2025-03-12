clc; clear;
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Parameters
scale_factor = 12.2; % pixels per cm
rod_length_cm = 0.7; % rod length in cm (7 mm)
rod_diameter_cm = 0.2; % rod diameter in cm (2 mm)
rod_length_px = rod_length_cm * scale_factor;
rod_diameter_px = rod_diameter_cm * scale_factor;
rod_half_width_px = (rod_diameter_cm / 2) * scale_factor;

rod_internal_resolution=100;

% Rotation parameters
rotation_rate_deg_per_s = 480/10; % degrees per second
rotation_rate_rad_per_s = rotation_rate_deg_per_s * (pi / 180); % rad/s
fps = 30; % frames per second
duration_s = 20; % total duration in seconds
t = linspace(0, duration_s, fps * duration_s); % time vector
angles = rotation_rate_rad_per_s * t; % Precomputed for efficiency % angles in radians over time

% Canvas size and rod position
canvas_size = [20, 20]; % [height, width] in pixels
center = canvas_size / 2;

% Generate synthetic images
num_frames = length(t);
images = uint8(zeros([canvas_size, num_frames]));


x_rot_all = zeros(num_frames, rod_internal_resolution);
y_rot_all = zeros(num_frames, rod_internal_resolution);

for i = 1:num_frames
    % Create a blank image
    img = zeros(canvas_size);

    % Define the rod as a line rotated by the current angle
    theta = angles(i); % current angle in radians
    x = linspace(-rod_length_px / 2, rod_length_px / 2, rod_internal_resolution); % rod length
    y = zeros(size(x));

    % Rotate the rod
    x_rot = cos(theta) * x - sin(theta) * y;
    y_rot = sin(theta) * x + cos(theta) * y;

    % Shift to center of canvas
    x_rot = x_rot + center(2);
    y_rot = y_rot + center(1);

    x_rot_all(i, :) = x_rot;
    y_rot_all(i, :) = y_rot;


    % Draw the rod on the image, activating pixels overlapping the rod region
    [x_pixel_grid, y_pixel_grid] = meshgrid(1:canvas_size(2), 1:canvas_size(1));

    % Calculate distances from each pixel to the line defined by the rod
    dx = x_rot(end) - x_rot(1);
    dy = y_rot(end) - y_rot(1);
    numerator = abs(dy .* (x_pixel_grid - x_rot(1)) - dx .* (y_pixel_grid - y_rot(1)));
    denominator = sqrt(dx^2 + dy^2);
    distance = numerator / denominator;

    % Find pixels that are within half the rod diameter and within rod endpoints
    dot_product1 = (x_pixel_grid - x_rot(1)) * dx + (y_pixel_grid - y_rot(1)) * dy;
    dot_product2 = (x_pixel_grid - x_rot(end)) * dx + (y_pixel_grid - y_rot(end)) * dy;
    within_endpoints = (dot_product1 >= 0) & (dot_product2 <= 0);

    % Activate pixels
    img((distance <= rod_half_width_px) & within_endpoints) = 1;


    images(:, :, i) = img;
end

% Extract orientation using regionprops
orientations = zeros(1, num_frames);
for i = 1:num_frames
    if any(images(:, :, i), 'all')
        props = regionprops(images(:, :, i), 'Orientation');
        orientations(i) = props.Orientation;
    else
        orientations(i) = NaN; % Handle cases where no objects are detected
    end
end


% correct orientations
orientations=correct_orientations(orientations);

% filter orientation
orientations = filter_orientations(orientations, 'movmean', 9);


% Compute rotation rate
rotation_rate = abs(compute_rotation_rate(t, orientations, 'central'));


%% Plot with subplots
figure(1); clf; set(gcf,'Position',[180 70 1000 700])
for i = 1:num_frames

    % Orientation over time
    subplot(2, 2, 1); cla;hold on;
    plot(t, orientations, 'k-o','LineWidth',1);
    plot(t(i), orientations(i), 'r.', 'MarkerSize', 20);
    grid on; box on; hold off;
    xlabel('Time (s)');
    ylabel('Orientation (deg)');
    title('Rod orientation from regionprops over time');

    % Rotation rate over time
    subplot(2, 2, 2); cla;hold on;
    plot(t, rotation_rate, 'k-o','LineWidth',1,'DisplayName', 'Measured rotation rate');
    plot(t(i), rotation_rate(i), 'r.', 'MarkerSize', 20, 'DisplayName', 'Measured current rotation rate');
    yline(rotation_rate_deg_per_s, 'b-', 'LineWidth', 2, 'DisplayName', 'Input rotation rate','Alpha',1);
    ylim([-5*rotation_rate_deg_per_s 5*rotation_rate_deg_per_s])
    grid on; box on; hold off;
    xlabel('Time (s)');
    ylabel('Absolute rotation rate (deg/s)');
    title('Rod rotation rate over time');
    legend('show', 'Location', 'best');
    %set(gca,'YScale','log')
    
    grid on;
    hold off;

    % Rotation rate error over time
    subplot(2, 2, 4); cla;hold on;
    plot(t, (rotation_rate-rotation_rate_deg_per_s)/rotation_rate_deg_per_s * 100, 'k-o','LineWidth',1,'DisplayName', 'Measured relative rotation rate error');
    plot(t(i), (rotation_rate(i)-rotation_rate_deg_per_s)/rotation_rate_deg_per_s * 100, 'r.', 'MarkerSize', 20, 'DisplayName', 'Measured current relative rotation rate error');
    mean_error = mean(abs((rotation_rate-rotation_rate_deg_per_s)/rotation_rate_deg_per_s * 100));
    yline(mean_error, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean relative error','Alpha',1);
    ylim([-2*mean_error 2* mean_error])
    grid on; box on; hold off;
    xlabel('Time (s)');
    ylabel('Relative rotation rate error ($\%$)');
    title('Relative rotation rate error over time');
    legend('show', 'Location', 'best');
    %set(gca,'YScale','log')
    
    grid on;
    hold off;

    % Rod image
    subplot(2, 2, 3); hold on;
    imagesc(images(:, :, i));

    % Draw pixel edges as vertical and horizontal lines
    for x_edge = 0.5:1:canvas_size(2)+0.5
        plot([x_edge, x_edge], [0.5, canvas_size(1)+0.5], 'k-', 'LineWidth', 0.5);
    end
    for y_edge = 0.5:1:canvas_size(1)+0.5
        plot([0.5, canvas_size(2)+0.5], [y_edge, y_edge], 'k-', 'LineWidth', 0.5);
    end


    % Draw the rod as a filled polygon representing its thickness
    dx = x_rot_all(i, end) - x_rot_all(i, 1);
    dy = y_rot_all(i, end) - y_rot_all(i, 1);
    length_factor = sqrt(dx^2 + dy^2);
    perp_x = -dy / length_factor * rod_half_width_px;
    perp_y = dx / length_factor * rod_half_width_px;

    % Define corners of the rectangle
    x_poly = [x_rot_all(i, 1) + perp_x, x_rot_all(i, end) + perp_x, x_rot_all(i, end) - perp_x, x_rot_all(i, 1) - perp_x];
    y_poly = [y_rot_all(i, 1) + perp_y, y_rot_all(i, end) + perp_y, y_rot_all(i, end) - perp_y, y_rot_all(i, 1) - perp_y];

    % Draw the filled polygon
    fill(x_poly, y_poly, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1);
    plot(x_rot_all(i,:),y_rot_all(i,:),'k--','LineWidth',2);

    box on;  daspect([1 1 1]); xlim([1,canvas_size(1)]); ylim([1,canvas_size(2)]);
    xticks(1:1:canvas_size(1)); yticks(1:1:canvas_size(2));
    title(sprintf('Frame %d / %d, Time: %.2f s', i, num_frames, t(i)));
    hold off
    pause(0.2);

end




%% functions
function rotation_rate = compute_rotation_rate(t, orientations, scheme)
% Computes the rotation rate using various finite difference schemes
% Inputs:
%   t - Time vector
%   orientations - Orientation angles
%   scheme - Finite difference scheme ('forward', 'backward', 'central', '5-point', '7-point')
% Output:
%   rotation_rate - Computed rotation rate

num_frames = length(t);
rotation_rate = zeros(1, num_frames);

switch scheme
    case 'forward'
        for i = 1:num_frames-1
            dt = t(i+1) - t(i);
            if ~isnan(orientations(i)) && ~isnan(orientations(i+1))
                rotation_rate(i) = (orientations(i+1) - orientations(i)) / dt;
            else
                rotation_rate(i) = NaN;
            end
        end
    case 'backward'
        for i = 2:num_frames
            dt = t(i) - t(i-1);
            if ~isnan(orientations(i)) && ~isnan(orientations(i-1))
                rotation_rate(i) = (orientations(i) - orientations(i-1)) / dt;
            else
                rotation_rate(i) = NaN;
            end
        end
    case 'central'
        for i = 2:num_frames-1
            dt = t(i+1) - t(i-1);
            if ~isnan(orientations(i-1)) && ~isnan(orientations(i+1))
                rotation_rate(i) = (orientations(i+1) - orientations(i-1)) / dt;
            else
                rotation_rate(i) = NaN;
            end
        end
    case '5-point'
        for i = 3:num_frames-2
            dt = t(i+2) - t(i-2);
            if ~any(isnan(orientations(i-2:i+2)))
                rotation_rate(i) = (-orientations(i+2) + 8*orientations(i+1) - 8*orientations(i-1) + orientations(i-2)) / (12*dt);
            else
                rotation_rate(i) = NaN;
            end
        end
    case '7-point'
        for i = 4:num_frames-3
            dt = t(i+3) - t(i-3);
            if ~any(isnan(orientations(i-3:i+3)))
                rotation_rate(i) = (orientations(i+3) - 9*orientations(i+2) + 45*orientations(i+1) - 45*orientations(i-1) + 9*orientations(i-2) - orientations(i-3)) / (60*dt);
            else
                rotation_rate(i) = NaN;
            end
        end
end
end

function corrected_orientations = correct_orientations(orientations)
% Corrects the orientations to remove jumps at +90 or -90 degrees and wraps properly
% Inputs:
%   orientations - Array of orientation angles over time
% Output:
%   corrected_orientations - Corrected orientation angles

corrected_orientations = orientations;
num_frames = length(orientations);

for frame = 2:num_frames
    % Skip if current or previous frame is NaN
    if isnan(corrected_orientations(frame)) || isnan(corrected_orientations(frame - 1))
        continue;
    end

    % Calculate the angle difference
    angleDifference = corrected_orientations(frame) - corrected_orientations(frame - 1);
    
    % Correct for wrapping beyond ±90 degrees
    while angleDifference > 90
        angleDifference = angleDifference - 180;
    end
    while angleDifference < -90
        angleDifference = angleDifference + 180;
    end
    
    % Update the corrected orientation
    corrected_orientations(frame) = corrected_orientations(frame - 1) + angleDifference;
end

end



function filtered_orientations = filter_orientations(orientations, method, window_size)
% Filters the orientations over time using MATLAB's smoothdata function
% Inputs:
%   orientations - Array of orientation angles over time
%   method - Smoothing method (e.g., 'movmean', 'gaussian', 'lowess')
%   window_size - Window size for the smoothing operation
% Output:
%   filtered_orientations - Filtered orientation angles

% Apply the smoothing filter
filtered_orientations = smoothdata(orientations, method,window_size);
end
