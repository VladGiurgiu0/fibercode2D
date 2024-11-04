% Define the path to the main folder
mainFolder = 'E:\ML_paper\Re1000_fiber_v4\Re1000_Fiber_3mm_01\';

% Get a list of all subfolders in the main folder
subfolders = dir(fullfile(mainFolder, 'loop=*'));

% Loop through each subfolder and rename it
for k = 1:length(subfolders)
    if subfolders(k).isdir
        % Get the current folder name
        oldFolderName = subfolders(k).name;
        
        % Extract the loop number from the folder name
        loopNum = sscanf(oldFolderName, 'loop=%d');
        
        % Calculate the new loop number
        newLoopNum = loopNum + 100;
        
        % Create the new folder name
        newFolderName = sprintf('loop=%d', newLoopNum);
        
        % Rename the folder
        movefile(fullfile(mainFolder, oldFolderName), fullfile(mainFolder, newFolderName));
    end
end
