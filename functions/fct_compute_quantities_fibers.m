function [] = fct_compute_quantities_fibers(p)
%% load tracked fibers
load([p.main_folder,'only_fibers/tracked_fibers.mat'])

%% remove trajectories which are shorter than given length
% based on the trajectory length
trajectory_length=sum(~isnan(positionsX),2);
positionsX(trajectory_length<p.min_track_length,:)=NaN;
positionsY(trajectory_length<p.min_track_length,:)=NaN;
orientationsMatrix(trajectory_length<p.min_track_length,:)=NaN;
majorAxisLengthsMatrix(trajectory_length < p.min_track_length, :) = NaN;
minorAxisLengthsMatrix(trajectory_length < p.min_track_length, :) = NaN;

%% process trajectories
positionsX_filtered = NaN(size(positionsX,1),size(positionsX,2));
positionsY_filtered = NaN(size(positionsX,1),size(positionsX,2));
orientationsMatrix_filtered = NaN(size(positionsX,1),size(positionsX,2));
velocityX = NaN(size(positionsX,1),size(positionsX,2));
velocityY = NaN(size(positionsX,1),size(positionsX,2));
omega_Z = NaN(size(positionsX,1),size(positionsX,2));

for i=1:size(positionsY,1)

    time = find(~isnan(orientationsMatrix(i, :)));
    if isempty(time)
    else
        % 1. filter
        [positionsX_filtered(i,time),~]=fit_data(time,positionsX(i,time)',p.filter_type,p.filter_kernel_positions);
        [positionsY_filtered(i,time),~]=fit_data(time,positionsY(i,time)',p.filter_type,p.filter_kernel_positions);
    
        [orientationsMatrix_filtered(i,time),~]=fit_data(time,orientationsMatrix(i,time)',p.filter_type,p.filter_kernel_orientations);
    
        % 2. compute derivative 
        % velocity of the center of mass in x
        velocityX(i,time)=compute_derivative(time,positionsX_filtered(i,time),p.derivation_scheme,0);
        % velocity of the center of mass in y
        velocityY(i,time)=compute_derivative(time,positionsY_filtered(i,time),p.derivation_scheme,0);
        % angular velocity around z - omega_fiber
        omega_Z(i,time)=compute_derivative(time,orientationsMatrix_filtered(i,time),p.derivation_scheme,0);
    end

    disp(['[1/1] Processed trajectory number: ',num2str(i),' out of ', num2str(size(positionsY,1))])

end

%% remove outlier trajectories
if p.remove_outliers==1
    % based on std of majoraxislength
    std_length = std(majorAxisLengthsMatrix,0,2,'omitnan');
    if p.plot_images==1
        figure();histogram(std_length,20); xlabel('std($L_1$)'); ylabel('counts'); 
        title('Histogram of standard deviation of the major axis length BEFORE REMOVAL')
    end
    positionsX_filtered(std_length>p.limit_std_length,:)=NaN;
    positionsY_filtered(std_length>p.limit_std_length,:)=NaN;
    orientationsMatrix_filtered(std_length>p.limit_std_length,:)=NaN;
    velocityX(std_length>p.limit_std_length,:)=NaN;
    velocityY(std_length>p.limit_std_length,:)=NaN;
    omega_Z(std_length>p.limit_std_length,:)=NaN;
    majorAxisLengthsMatrix(std_length>p.limit_std_length,:)=NaN;
    minorAxisLengthsMatrix(std_length>p.limit_std_length,:)=NaN;
    
    remaining_trajectories = mean(velocityX,2,'omitnan');
    
    % based on the size of the minoraxislength
    mean_diameter = mean(minorAxisLengthsMatrix,2,'omitnan');
    if p.plot_images==1
        figure();histogram(mean_diameter,20); xlabel('$L_2$'); ylabel('counts'); 
        title('Histogram of minor axis length BEFORE REMOVAL')
    end
    positionsX_filtered(mean_diameter>p.max_fiber_diameter,:)=NaN;
    positionsY_filtered(mean_diameter>p.max_fiber_diameter,:)=NaN;
    orientationsMatrix_filtered(mean_diameter>p.max_fiber_diameter,:)=NaN;
    velocityX(mean_diameter>p.max_fiber_diameter,:)=NaN;
    velocityY(mean_diameter>p.max_fiber_diameter,:)=NaN;
    omega_Z(mean_diameter>p.max_fiber_diameter,:)=NaN;
    majorAxisLengthsMatrix(mean_diameter>p.max_fiber_diameter,:)=NaN;
    minorAxisLengthsMatrix(mean_diameter>p.max_fiber_diameter,:)=NaN;
    
    % based on std of velocity
    std_velocityX = std(velocityX,0,2,'omitnan');
    if p.plot_images==1
        figure();histogram(std_velocityX,20); xlabel('std($u_x$)'); ylabel('counts'); 
        title('Histogram of standard deviation of stream-wise velocity within each trajectory BEFORE REMOVAL')
    end
    positionsX_filtered(std_velocityX>p.limit_std_velocityX,:)=NaN;
    positionsY_filtered(std_velocityX>p.limit_std_velocityX,:)=NaN;
    orientationsMatrix_filtered(std_velocityX>p.limit_std_velocityX,:)=NaN;
    velocityX(std_velocityX>p.limit_std_velocityX,:)=NaN;
    velocityY(std_velocityX>p.limit_std_velocityX,:)=NaN;
    omega_Z(std_velocityX>p.limit_std_velocityX,:)=NaN;
    majorAxisLengthsMatrix(std_velocityX>p.limit_std_velocityX,:)=NaN;
    minorAxisLengthsMatrix(std_velocityX>p.limit_std_velocityX,:)=NaN;
    
    
    % based on wall-normal location
    % remove everything below the wall
    % positionsX_filtered(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % positionsY_filtered(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % orientationsMatrix_filtered(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % velocityX(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % velocityY(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % omega_Z(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % majorAxisLengthsMatrix(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
    % minorAxisLengthsMatrix(mean(positionsY,2,'omitnan')>wall_position,:)=NaN;
end

%% remove trajectories which don't contain data
does_not_contain_data = isnan(mean(velocityX,2,"omitnan"));

positionsX_filtered(does_not_contain_data,:)=[];
positionsY_filtered(does_not_contain_data,:)=[];
orientationsMatrix_filtered(does_not_contain_data,:)=[];
velocityX(does_not_contain_data,:)=[];
velocityY(does_not_contain_data,:)=[];
omega_Z(does_not_contain_data,:)=[];
majorAxisLengthsMatrix(does_not_contain_data, :) = [];
minorAxisLengthsMatrix(does_not_contain_data, :) = [];

%% save trajectory data
save([p.main_folder,'only_fibers/','quantities_fibers.mat'], ...
    'positionsX_filtered','positionsY_filtered','velocityX',"velocityY","omega_Z",'orientationsMatrix_filtered',...
    'majorAxisLengthsMatrix','minorAxisLengthsMatrix','p','-v7.3');