clc; clear; close all;
% Define the example to use (1 to 6)
example = 5;

% Define grid
x = -10:10;
y = -10:10;
[X, Y] = meshgrid(x, y);

% Initialize U and V
U = zeros(size(X));
V = zeros(size(Y));

% Analytical flow field parameters
Omega = 1;  % for solid body rotation
Gamma = 1;  % for shear flow
U0 = 1;     % for uniform flow
V0 = 0;     % for uniform flow
Q = 1;      % for source flow
alpha = 1;  % for hyperbolic flow

% Define analytical flow fields
switch example
    case 1
        % Solid Body Rotation
        U = -Omega * Y;
        V = Omega * X;
        analytical_vorticity = 2 * Omega * ones(size(U));
        analytical_lambda_ci = abs(Omega) * ones(size(U));
        analytical_strain_rate = zeros(size(U));
    case 2
        % Simple Shear Flow
        U = Gamma * Y;
        V = zeros(size(Y));
        analytical_vorticity = Gamma * ones(size(U));
        analytical_lambda_ci = zeros(size(U));
        analytical_strain_rate = (Gamma / 2) * ones(size(U));
    case 3
        % Uniform Flow
        U = U0 * ones(size(X));
        V = V0 * ones(size(Y));
        analytical_vorticity = zeros(size(U));
        analytical_lambda_ci = zeros(size(U));
        analytical_strain_rate = zeros(size(U));
    case 4
        % Source Flow
        U = (Q / (2 * pi)) * (X ./ (X.^2 + Y.^2));
        V = (Q / (2 * pi)) * (Y ./ (X.^2 + Y.^2));
        analytical_vorticity = zeros(size(U));
        analytical_lambda_ci = zeros(size(U));
        analytical_strain_rate = Q ./ (2 * pi * sqrt(X.^2 + Y.^2));
    case 5
        % Vortex Flow
        U = - (Gamma / (2 * pi)) * (Y ./ (X.^2 + Y.^2));
        V = (Gamma / (2 * pi)) * (X ./ (X.^2 + Y.^2));
        analytical_vorticity = (Gamma ./ (pi * (X.^2 + Y.^2)));
        analytical_lambda_ci = (Gamma ./ (2 * pi * (X.^2 + Y.^2)));
        analytical_strain_rate = zeros(size(U));
    case 6
        % Hyperbolic Flow
        U = alpha * X;
        V = -alpha * Y;
        analytical_vorticity = zeros(size(U));
        analytical_lambda_ci = zeros(size(U));
        analytical_strain_rate = alpha * ones(size(U));
    otherwise
        error('Invalid example number. Choose between 1 and 6.');
end

% Compute the velocity gradient tensor
gradTensor = compute_velocity_gradient_tensor(U, V);

% Compute vorticity
vorticity = compute_vorticity(gradTensor);

% Compute swirling strength
lambda_ci = compute_swirling_strength(gradTensor);

% Compute strain rate
strain_rate = compute_strain_rate(gradTensor);

% Plot results
figure;

subplot(3, 3, 1);
imagesc(x, y, U);
colorbar;
title('U velocity');

subplot(3, 3, 2);
imagesc(x, y, V);
colorbar;
title('V velocity');

subplot(3, 3, 3);
quiver(X, Y, U, V);
title('Velocity Field');

subplot(3, 3, 4);
imagesc(x, y, vorticity);
colorbar;
title('Computed Vorticity');

subplot(3, 3, 5);
imagesc(x, y, lambda_ci);
colorbar;
title('Computed Swirling Strength');

subplot(3, 3, 6);
imagesc(x, y, strain_rate);
colorbar;
title('Computed Strain Rate');

subplot(3, 3, 7);
imagesc(x, y, analytical_vorticity);
colorbar;
title('Analytical Vorticity');

subplot(3, 3, 8);
imagesc(x, y, analytical_lambda_ci);
colorbar;
title('Analytical Swirling Strength');

subplot(3, 3, 9);
imagesc(x, y, analytical_strain_rate);
colorbar;
title('Analytical Strain Rate');

% Display differences
figure;

subplot(1, 3, 1);
imagesc(x, y, vorticity - analytical_vorticity);
colorbar;
title('Vorticity Error');

subplot(1, 3, 2);
imagesc(x, y, lambda_ci - analytical_lambda_ci);
colorbar;
title('Swirling Strength Error');

subplot(1, 3, 3);
imagesc(x, y, strain_rate - analytical_strain_rate);
colorbar;
title('Strain Rate Error');
