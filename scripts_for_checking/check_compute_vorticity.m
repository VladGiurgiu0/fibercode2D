% Test the computeVorticity function with a known velocity field
clc; clear; close all;

%% check with solid body rotation
clc; clear; close all;
% Define the grid
x = 1:100;
y = 1:50;
[xx, yy] = meshgrid(x, y);

% Define the velocity field
% Example: Solid body rotation
omega = 2; % angular velocity
u = -omega * yy;
v = omega * xx;

% Define the expected vorticity for solid body rotation
expectedVorticity = 2 * omega * ones(size(xx));

% Compute vorticity using the function
%computedVorticity = computeVorticity(u, v);

gradTensor = compute_velocity_gradient_tensor(u, v);
computedVorticity = compute_vorticity(gradTensor);

% Plot the results
figure;
subplot(1, 3, 1);
quiver(xx,yy,u,v,1);
axis equal
%colorbar;
title('Expected Vorticity');
xlabel('X');
ylabel('Y');

subplot(1, 3, 2);
imagesc(expectedVorticity);
colorbar;
axis equal
title('Expected Vorticity');
xlabel('X');
ylabel('Y');

subplot(1, 3, 3);
imagesc(computedVorticity);
colorbar;
axis equal
title('Computed Vorticity');
xlabel('X');
ylabel('Y');

% Calculate and display the error
error = norm(computedVorticity - expectedVorticity, 'fro') / norm(expectedVorticity, 'fro');
fprintf('Relative Frobenius norm error: %e\n', error);

%% check with analytical solution
clc; clear; close all;
% Define the grid
x = 1:100;
y = 1:50;
[xx, yy] = meshgrid(x, y);

u = xx.^2 + yy.^2;
v = xx.*yy + xx.^3;

expectedVorticity = -yy + 2*xx.^2;

% Compute vorticity using the function
%computedVorticity = computeVorticity(u, v);


gradTensor = compute_velocity_gradient_tensor(u, v);
computedVorticity = compute_vorticity(gradTensor);

% Plot the results
figure;
subplot(1, 3, 1);
quiver(xx,yy,u,v,1);
axis equal
%colorbar;
title('Expected Vorticity');
xlabel('X');
ylabel('Y');


subplot(1, 3, 2);
imagesc(expectedVorticity);
colorbar;
axis equal
title('Expected Vorticity');
xlabel('X');
ylabel('Y');
clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

subplot(1, 3, 3);
imagesc(computedVorticity);
colorbar;
axis equal
title('Computed Vorticity');
xlabel('X');
ylabel('Y');
clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

% Calculate and display the error
error = norm(computedVorticity - expectedVorticity, 'fro') / norm(expectedVorticity, 'fro');
fprintf('Relative Frobenius norm error: %e\n', error);

%% check with analytical solution 2
clc; clear; close all;
% Define the grid
x = 1:100;
y = 1:50;
[xx, yy] = meshgrid(x, y);

u = xx + yy.^2;
v = 2*xx + 2*yy;

expectedVorticity = 2 - 2*yy;

% Compute vorticity using the function
%computedVorticity = computeVorticity(u, v);

gradTensor = compute_velocity_gradient_tensor(u, v);
computedVorticity = compute_vorticity(gradTensor);

% Plot the results
figure;
subplot(1, 3, 1);
quiver(xx,yy,u,v,1);
axis equal
%colorbar;
title('Expected Vorticity');
xlabel('X');
ylabel('Y');

subplot(1, 3, 2);
imagesc(expectedVorticity);
colorbar;
axis equal
title('Expected Vorticity');
xlabel('X');
ylabel('Y');
clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

subplot(1, 3, 3);
imagesc(computedVorticity);
colorbar;
axis equal
title('Computed Vorticity');
xlabel('X');
ylabel('Y');
clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

% Calculate and display the error
error = norm(computedVorticity - expectedVorticity, 'fro') / norm(expectedVorticity, 'fro');
fprintf('Relative Frobenius norm error: %e\n', error);

%% check with analytical solution 3
clc; clear; close all;
% Define the grid
x = 1:100;
y = 1:50;
[xx, yy] = meshgrid(x, y);

u = 2*yy;
v = 3*xx;

expectedVorticity = 1;

% Compute vorticity using the function
%computedVorticity = computeVorticity(u, v);

gradTensor = compute_velocity_gradient_tensor(u, v);
computedVorticity = compute_vorticity(gradTensor);

% Plot the results
figure;
subplot(1, 3, 1);
quiver(xx,yy,u,v,1);
axis equal
%colorbar;
title('Expected Vorticity');
xlabel('X');
ylabel('Y');


subplot(1, 3, 2);
imagesc(expectedVorticity);
colorbar;
axis equal
title('Expected Vorticity');
xlabel('X');
ylabel('Y');
%clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

subplot(1, 3, 3);
imagesc(computedVorticity);
colorbar;
axis equal
title('Computed Vorticity');
xlabel('X');
ylabel('Y');
%clim([min(expectedVorticity(:)),max(expectedVorticity(:))])

% Calculate and display the error
error = norm(computedVorticity - expectedVorticity, 'fro') / norm(expectedVorticity, 'fro');
fprintf('Relative Frobenius norm error: %e\n', error);



%% check with PIV data
load('/Users/vlad/Owncloud/Research/Papers/machine learning/data/Re1000_fiber/only_tracers/output_PaIRS/out_0000.mat')

vorticity = computeVorticity(U,V);

figure;
subplot(1,3,1);
imagesc(U);
daspect([1 1 1]); colorbar; box on;

subplot(1,3,2);
imagesc(V);
daspect([1 1 1]); colorbar; box on;

subplot(1,3,3);
imagesc(vorticity);
daspect([1 1 1]); colorbar; box on;