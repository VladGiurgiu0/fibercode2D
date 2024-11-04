clc; clear; close all

% Base directories
source_base_dir = 'D:\ML_paper\Re1000_fiber_v4\Re1000_Fiber_3mm\loop=';
destination_base_dir = 'D:\ML_paper\Re1000_fiber_v4\Re1000_Fiber_3mm_only_images\loop=';

% Loop through each loop directory
for loop_num = 0:198
    % Construct the full source and destination paths
    source_dir = [source_base_dir, num2str(loop_num), '\only_tracers\'];
    destination_dir = [destination_base_dir, num2str(loop_num), '\only_tracers\'];

    % Create the destination directory if it does not exist
    if ~exist(destination_dir, 'dir')
        mkdir(destination_dir);
    end

    % Loop through each image from tracers_0001.tif to tracers_0300.tif
    for img_num = 1:300
        % Create the image file name with leading zeros
        img_name = sprintf('tracers_%04d.tif', img_num);

        % Full path to the source and destination images
        source_file = fullfile(source_dir, img_name);
        destination_file = fullfile(destination_dir, img_name);

        % Copy the file from source to destination
        if exist(source_file, 'file')
            copyfile(source_file, destination_file);
        else
            fprintf('Warning: %s does not exist.\n', source_file);
        end
    end
end

disp('Image copying process completed.');
