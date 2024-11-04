clear; close all; clc;
set(0,'DefaultAxesFontSize',20,'DefaultFigureColor', [1 1 1],'defaultfigureposition',[50 100 1800 600])
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Input
fps = 1478;             % aquisition frequnecy [Hz]
px_mm = 35.36;          % scale factor [px/mm]
search_radius = 50;     % radius within which to consider the flow field [px]

nr_loops = 1;
nr_frames = 5;
root_folder = 'example_data/';

cont_fiber = 0;
rotation_pdf = [];
%%
for iiii=1:nr_loops
% Load data
disp(['Loop number computing: ', num2str(iiii)])
main_folder = [root_folder,'loop=',num2str(iiii-1,'%d'),'/'];
% Fibers:
load([main_folder,'/only_fibers/','mask_fibers.mat']);
load([main_folder,'/only_fibers/','quantities_fibers.mat'])
%load([main_folder,'/only_fibers/','tracked_fibers.mat'])
% Flow:
for jjjj=1:nr_frames
    load([main_folder,'/only_tracers/out_PaIRS/','flow_',num2str(jjjj-1,'%04d'),'.mat']);
    U_IMG_temp(:,:,jjjj) = U_IMG;
    V_IMG_temp(:,:,jjjj) = V_IMG;
    lambda_ci_IMG_temp(:,:,jjjj) = lambda_ci_IMG;
    vorticity_IMG_temp(:,:,jjjj) = vorticity_IMG;
    shear_rate_IMG_temp(:,:,jjjj) = shear_rate_IMG;

    %disp(['Loaded ',num2str(jjjj)])
end

%% Gather statistics
[x, y] = meshgrid(linspace(1,size(U_IMG_temp,2),size(U_IMG_temp,2)),linspace(1,size(U_IMG_temp,1),size(U_IMG_temp,1)));
U_field = U_IMG_temp./px_mm.*fps./1000; % (m/s)
V_field = V_IMG_temp./px_mm.*fps./1000; % (m/s)

lambda_field = lambda_ci_IMG_temp.*fps^2;
vorticity_field = vorticity_IMG_temp.*fps;
shear_field = shear_rate_IMG_temp.*fps;

final_fiber_cont = cont_fiber;

for fiber_Id = 1+cont_fiber:size(omega_Z,1)+cont_fiber
    valid_frames = find(~isnan(omega_Z(fiber_Id-final_fiber_cont,1:end-1)));
    cont = 1;   
    for i = valid_frames
        center = round([positionsY_filtered(fiber_Id-final_fiber_cont,i), positionsX_filtered(fiber_Id-final_fiber_cont,i)]);
        x_range = max(1, center(2)-search_radius):min(size(U_IMG_temp,2), center(2)+search_radius);
        y_range = max(1, center(1)-search_radius):min(size(U_IMG_temp,1), center(1)+search_radius);
        [Y, X] = meshgrid(y_range, x_range);
        inside_circle = ((X - center(2)).^2 + (Y - center(1)).^2) <= search_radius^2;
        inside_circle = inside_circle';
        masked_flow_field = imcomplement(fiberMask_cell{i});
        masked_flow_field = double(masked_flow_field(y_range, x_range));   
        indices = (masked_flow_field == 1);
        masked_flow_field(indices) = 0;        
        masked_flow_field(~indices) = NaN;

        % y [m]
        local(fiber_Id).y(cont) = center(1)/(px_mm*1000);
        % U-velocity [m/s]
        coupled_flow_field = U_field(y_range, x_range, i) + masked_flow_field;
        coupled_flow_field(~inside_circle) = NaN;
        local(fiber_Id).U(cont) = mean(coupled_flow_field,'all','OmitNaN');
        % V-velocity [m/s]
        coupled_flow_field = V_field(y_range, x_range, i) + masked_flow_field;
        coupled_flow_field(~inside_circle) = NaN;
        local(fiber_Id).V(cont) = mean(coupled_flow_field,'all','OmitNaN');
        % Lambda_ci [1/s^2]
        coupled_flow_field = lambda_field(y_range, x_range, i) + masked_flow_field;
        coupled_flow_field(~inside_circle) = NaN;
        local(fiber_Id).lambda(cont) = mean(coupled_flow_field,'all','OmitNaN');
        % Vorticity [1/s]
        coupled_flow_field = vorticity_field(y_range, x_range, i) + masked_flow_field;
        coupled_flow_field(~inside_circle) = NaN;
        local(fiber_Id).vorticity(cont) = mean(coupled_flow_field,'all','OmitNaN');
        % shear rate [1/s]
        coupled_flow_field = shear_field(y_range, x_range, i) + masked_flow_field;
        coupled_flow_field(~inside_circle) = NaN;
        local(fiber_Id).shear(cont) = mean(coupled_flow_field,'all','OmitNaN');
        % Rotation rate [rad/s]
        local(fiber_Id).rotation(cont) = deg2rad(omega_Z(fiber_Id-final_fiber_cont,i)).*fps;
        % Velocity [m/s]
        local(fiber_Id).u_fiber(cont) = velocityX(fiber_Id-final_fiber_cont,i).*(1/(px_mm*1000)).*fps;
        local(fiber_Id).v_fiber(cont) = velocityY(fiber_Id-final_fiber_cont,i).*(1/(px_mm*1000)).*fps;

        % Assign major and minor axis lengths [mm]
        local(fiber_Id).major_axis_length(cont) = majorAxisLengthsMatrix(fiber_Id-final_fiber_cont, i) / px_mm;
        local(fiber_Id).minor_axis_length(cont) = minorAxisLengthsMatrix(fiber_Id-final_fiber_cont, i) / px_mm;



        %%% flow quantities interpolated at the center of mass of the fiber
        % Interpolate U-velocity at the fiber's center of mass
        U_interpolated = interp2(U_field(:,:,i), center(2), center(1));
        local(fiber_Id).U_interpolated(cont) = U_interpolated;
        
        % Interpolate V-velocity at the fiber's center of mass
        V_interpolated = interp2(V_field(:,:,i), center(2), center(1));
        local(fiber_Id).V_interpolated(cont) = V_interpolated;
        
        % Interpolate Vorticity at the fiber's center of mass
        vorticity_interpolated = interp2(vorticity_field(:,:,i), center(2), center(1));
        local(fiber_Id).vorticity_interpolated(cont) = vorticity_interpolated;
        
        % Interpolate Shear rate at the fiber's center of mass
        shear_interpolated = interp2(shear_field(:,:,i), center(2), center(1));
        local(fiber_Id).shear_interpolated(cont) = shear_interpolated;
        
        % Interpolate Lambda_ci at the fiber's center of mass
        lambda_interpolated = interp2(lambda_field(:,:,i), center(2), center(1));
        local(fiber_Id).lambda_interpolated(cont) = lambda_interpolated;

        cont = cont + 1;
    end
    cont_fiber = cont_fiber + 1;
end
end

%% put data into matrixes and save
% Determine the maximum number of timesteps
max_timesteps = max(arrayfun(@(x) length(x.U), local));

% Initialize empty matrices for each field
U_data = NaN(length(local), max_timesteps);
V_data = NaN(length(local), max_timesteps);
lambda_data = NaN(length(local), max_timesteps);
vorticity_data = NaN(length(local), max_timesteps);
shear_data = NaN(length(local), max_timesteps);
rotation_data = NaN(length(local), max_timesteps);
fiber_y_data = NaN(length(local), max_timesteps);
fiber_U_data = NaN(length(local), max_timesteps);
fiber_V_data = NaN(length(local), max_timesteps);

major_axis_length_data = NaN(length(local), max_timesteps);
minor_axis_length_data = NaN(length(local), max_timesteps);

% Initialize empty matrices for interpolated data
U_interpolated_data = NaN(length(local), max_timesteps);
V_interpolated_data = NaN(length(local), max_timesteps);
vorticity_interpolated_data = NaN(length(local), max_timesteps);
shear_interpolated_data = NaN(length(local), max_timesteps);
lambda_interpolated_data = NaN(length(local), max_timesteps);

% Loop through each fiber in the local structure
for fiber_Id = 1:length(local)
    U_data(fiber_Id, 1:length(local(fiber_Id).U)) = local(fiber_Id).U;
    V_data(fiber_Id, 1:length(local(fiber_Id).V)) = local(fiber_Id).V;
    lambda_data(fiber_Id, 1:length(local(fiber_Id).lambda)) = local(fiber_Id).lambda;
    vorticity_data(fiber_Id, 1:length(local(fiber_Id).vorticity)) = local(fiber_Id).vorticity;
    shear_data(fiber_Id, 1:length(local(fiber_Id).shear)) = local(fiber_Id).shear;
    rotation_data(fiber_Id, 1:length(local(fiber_Id).rotation)) = local(fiber_Id).rotation;
    fiber_y_data(fiber_Id, 1:length(local(fiber_Id).y)) = local(fiber_Id).y;
    fiber_U_data(fiber_Id, 1:length(local(fiber_Id).u_fiber)) = local(fiber_Id).u_fiber;
    fiber_V_data(fiber_Id, 1:length(local(fiber_Id).v_fiber)) = local(fiber_Id).v_fiber;

    major_axis_length_data(fiber_Id, 1:length(local(fiber_Id).major_axis_length)) = local(fiber_Id).major_axis_length;
    minor_axis_length_data(fiber_Id, 1:length(local(fiber_Id).minor_axis_length)) = local(fiber_Id).minor_axis_length;


    % Assign interpolated data to matrices
    U_interpolated_data(fiber_Id, 1:length(local(fiber_Id).U_interpolated)) = local(fiber_Id).U_interpolated;
    V_interpolated_data(fiber_Id, 1:length(local(fiber_Id).V_interpolated)) = local(fiber_Id).V_interpolated;
    vorticity_interpolated_data(fiber_Id, 1:length(local(fiber_Id).vorticity_interpolated)) = local(fiber_Id).vorticity_interpolated;
    shear_interpolated_data(fiber_Id, 1:length(local(fiber_Id).shear_interpolated)) = local(fiber_Id).shear_interpolated;
    lambda_interpolated_data(fiber_Id, 1:length(local(fiber_Id).lambda_interpolated)) = local(fiber_Id).lambda_interpolated;
    
end

% Now U_data, V_data, lambda_data, vorticity_data, shear_data, rotation_data,
% and y_plus_data are matrices with rows corresponding to fiber_Id and columns to timesteps.
save([root_folder,'fibers_coupled_data_50px.mat'],...
    'U_data', 'V_data', 'lambda_data', 'vorticity_data', 'shear_data', ...
    'rotation_data', 'fiber_y_data', 'fiber_U_data', 'fiber_V_data', ...
    'major_axis_length_data', 'minor_axis_length_data', ...
    'U_interpolated_data', 'V_interpolated_data', 'vorticity_interpolated_data', ...
    'shear_interpolated_data', 'lambda_interpolated_data');

