clc; clear; close all
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath("functions/")

%% parameters
p.root_folder = 'example_data/';
p.nr_loops = 1;
p.plot_images = 0; % set to 1 to visually check results

% discretisation
p.threshold_for_binarization_keep_fibers=0.2; % 0 - 1  - this defines which fibers to keep
p.threshold_for_binarization_remove_fibers=0.05; % 0 - 1  - this defines which fibers to delete from the tracer images

p.length_threshold_keep_fibers = 30; % Fibers with major axis length above this are kept in the fiber images[px]
p.length_threshold_remove_fibers = 15; % Fibers with major axis length above this are deleted in the tracer images [px]
p.dilation_size=1; % how much to dilate fibers before removing them from tracer images [px]

% tracking
p.distanceThreshold = 30; % search radius [px]

% computation of fiber quantities
p.filter_type = 'rlowess'; % which filter type should be used to filter the x and y positions and orientations
p.filter_kernel_positions = 50; % filter kernel for the positions in time-steps [-]
p.filter_kernel_orientations = 50; % filter kernel for the positions in time-steps [-]
p.derivation_scheme = '2 points stencil'; % scheme for computing derivatives; options: '2 points stencil', '5 points stencil'
p.min_track_length = 3; % remove fiber trajectories shorter than this value in timesteps [-]
p.remove_outliers = 0; % 0 or 1 - if 1 then remove outliers detected based on the thresholds below
%wall_position = 0; % remove fiber trajectories with the mean wall-normal location above this value [px] 
p.limit_std_velocityX = 10; % remove fiber trajectories with a standard deviation of the stream-wise velocity within a trajectory above this value
p.limit_std_length = 30; % remove fiber trajectories with a standard deviation of the majoraxislength above this value
p.max_fiber_diameter = 15; % remove fiber trajectories with minoraxislength larger than this value

% computation of flow quantities
p.nr_frames=5; % how many frames should be computed
p.which_flow_field = 'pairs'; % 'pairs', 'raft'
% time filter
p.filter_in_time = 1;
p.filter_time_kernel=20;
p.filter_time_type = 'rlowess';
% if 'pairs'
p.interpolate_velocity = 1; % if 1 interpolate velocity field onto a 1x1 px grid and then compute derivatives
                          % if 0 then the velocity field and its derivatives are upsampled
                          % without interpolation -> window size is needed
p.interp_type = 'linear'; % 'linear','cubic', 'makima', 'spline'
p.IW_size = 16;
% spatial filtering
p.filter_velocity_field = 1; % if 1 then use a gaussian filter on the velocity field before computing derivatives
p.filter_sigma = 10/16; % sigma; kernel size is: 2*ceil(2*sigma)+1

%% compute
for iiii = 1:p.nr_loops
    p.main_folder = [p.root_folder,'loop=',num2str(iiii-1,'%d'),'/'];

    % save parameters
    save([p.main_folder,'parameters.mat'],"p");

    % discrimination
    fct_discriminate_fibers_from_tracers(p);

    % tracking
     fct_track_fibers(p);

    % compute fiber quantities
    fct_compute_quantities_fibers(p);

    disp(['----------------------------- Finished loop ',num2str(iiii)])

end

%% AFTER PIV: compute flow and quantities
for iiii = 1:p.nr_loops
    p.main_folder = [p.root_folder,'loop=',num2str(iiii-1,'%d'),'/'];
    fct_compute_flow(p,iiii);
    disp(['Finished loop nr: ',num2str(iiii)])
end
