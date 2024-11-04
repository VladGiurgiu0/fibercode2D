function [] = fct_track_fibers(p)
%% load images
% image data store
imds = imageDatastore([p.main_folder,'only_fibers/'],'FileExtensions','.tif');


% Load the sequence of preprocessed images (fibers only)
numFrames = size(imds.Files,1);

% Initialize arrays to store centroids, orientations, and trajectories
centroids = cell(numFrames, 1);
orientations = cell(numFrames, 1);
trajectories = cell(numFrames, 1);
pixelLists = cell(numFrames, 1); % Initialize to store pixel lists
majorAxisLengths = cell(numFrames, 1); % Initialize to store MajorAxisLength
minorAxisLengths = cell(numFrames, 1); % Initialize to store MinorAxisLength


% Process each frame to compute centroids and orientations
for frame = 1:numFrames
    binaryImg = imbinarize(readimage(imds,frame),0);

    % Measure properties of image regions
    stats = regionprops(binaryImg, 'Centroid', 'Orientation','PixelList', 'MajorAxisLength', 'MinorAxisLength');

    % Store centroids and orientations
    centroids{frame} = cat(1, stats.Centroid);
    orientations{frame} = cat(1, stats.Orientation);
    pixelLists{frame} = {stats.PixelList}; % Cell array to store PixelList
    majorAxisLengths{frame} = cat(1, stats.MajorAxisLength);
    minorAxisLengths{frame} = cat(1, stats.MinorAxisLength);

    disp(['[1/2] Processed image number: ',num2str(frame),' out of ', num2str(numFrames)])
end

% Initialize trajectories with the first frame
trackedCentroids = cell(numFrames, 1);
trackedOrientations = cell(numFrames, 1);
fiberIDs = (1:size(centroids{1}, 1))';
trackedCentroids{1} = [fiberIDs centroids{1}];
trackedOrientations{1} = [fiberIDs orientations{1}];
trackedPixelLists{1} = pixelLists{1};
trackedMajorAxisLengths = cell(numFrames, 1); % Initialize to track MajorAxisLength
trackedMinorAxisLengths = cell(numFrames, 1); % Initialize to track MinorAxisLength
trackedMajorAxisLengths{1} = [fiberIDs majorAxisLengths{1}];
trackedMinorAxisLengths{1} = [fiberIDs minorAxisLengths{1}];

% Store initial positions in trajectories
for i = 1:length(fiberIDs)
    trajectories{fiberIDs(i)} = centroids{1}(i, :);
end

% Nearest neighbor tracking across frames
for frame = 2:numFrames
    prevFrameCentroids = trackedCentroids{frame-1};
    currCentroids = centroids{frame};
    currOrientations = orientations{frame};
    currPixelLists = pixelLists{frame};
    currMajorAxisLengths = majorAxisLengths{frame};
    currMinorAxisLengths = minorAxisLengths{frame};
    
    % Initialize arrays for tracked centroids and orientations
    trackedCentroids{frame} = [];
    trackedOrientations{frame} = [];
    trackedPixelLists{frame} = cell(size(currCentroids, 1), 1);
    trackedMajorAxisLengths{frame} = [];
    trackedMinorAxisLengths{frame} = [];
    
    % Create a copy of fiberIDs to keep track of used IDs
    availableFiberIDs = fiberIDs;
    
    % Create a list to keep track of which current centroids are already matched
    matchedCurrIndices = false(size(currCentroids, 1), 1);
    
    % Track each fiber from the previous frame
    for i = 1:size(prevFrameCentroids, 1)
        prevID = prevFrameCentroids(i, 1);
        prevCentroid = prevFrameCentroids(i, 2:3);
        
        distances = sqrt(sum((currCentroids - prevCentroid).^2, 2));
        [minDistance, minIdx] = min(distances);
        
        if minDistance < p.distanceThreshold && ~matchedCurrIndices(minIdx)
            % Append the current centroid and orientation with the same ID
            trackedCentroids{frame} = [trackedCentroids{frame}; prevID, currCentroids(minIdx, :)];
            trackedOrientations{frame} = [trackedOrientations{frame}; prevID, currOrientations(minIdx)];
            trackedPixelLists{frame}{prevID} = currPixelLists{minIdx};
            trackedMajorAxisLengths{frame} = [trackedMajorAxisLengths{frame}; prevID, currMajorAxisLengths(minIdx)];
            trackedMinorAxisLengths{frame} = [trackedMinorAxisLengths{frame}; prevID, currMinorAxisLengths(minIdx)];
            
            % Store position in trajectories
            trajectories{prevID} = [trajectories{prevID}; currCentroids(minIdx, :)];
            
            % Mark the current centroid as matched
            matchedCurrIndices(minIdx) = true;
            
            % Remove the matched fiber from availableFiberIDs
            availableFiberIDs(availableFiberIDs == prevID) = [];
        end
    end

    % Assign new IDs to new fibers that appear in the current frame
    unmatchedCurrIndices = find(~matchedCurrIndices);
    for k = 1:length(unmatchedCurrIndices)
        newFiberIdx = unmatchedCurrIndices(k);
        newID = max(fiberIDs) + 1;
        fiberIDs = [fiberIDs; newID];
        trackedCentroids{frame} = [trackedCentroids{frame}; newID, currCentroids(newFiberIdx, :)];
        trackedOrientations{frame} = [trackedOrientations{frame}; newID, currOrientations(newFiberIdx)];
        trackedPixelLists{frame}{newID} = currPixelLists{newFiberIdx};
        trackedMajorAxisLengths{frame} = [trackedMajorAxisLengths{frame}; newID, currMajorAxisLengths(newFiberIdx)];
        trackedMinorAxisLengths{frame} = [trackedMinorAxisLengths{frame}; newID, currMinorAxisLengths(newFiberIdx)];
        
        
        % Initialize trajectory for the new fiber
        trajectories{newID} = currCentroids(newFiberIdx, :);
    end
    disp(['[2/2] Tracking on image number: ',num2str(frame),' out of ', num2str(numFrames)])
end


%% reshape the data
% Initialize matrices to store positions and orientations
numFibers = max(fiberIDs);
positionsX = NaN(numFibers, numFrames);
positionsY = NaN(numFibers, numFrames);
orientationsMatrix = NaN(numFibers, numFrames);
pixelListsMatrix = cell(numFibers, numFrames);
majorAxisLengthsMatrix = NaN(numFibers, numFrames);
minorAxisLengthsMatrix = NaN(numFibers, numFrames);

% Populate the matrices
for frame = 1:numFrames
    for i = 1:size(trackedCentroids{frame}, 1)
        fiberID = trackedCentroids{frame}(i, 1);
        centroid = trackedCentroids{frame}(i, 2:3);
        orientation = trackedOrientations{frame}(i, 2);
        pixelList = trackedPixelLists{frame}{fiberID};
        majorAxisLength = trackedMajorAxisLengths{frame}(i, 2);
        minorAxisLength = trackedMinorAxisLengths{frame}(i, 2);
        
        % Store the positions and orientations in the matrices
        positionsX(fiberID, frame) = centroid(1);
        positionsY(fiberID, frame) = centroid(2);
        orientationsMatrix(fiberID, frame) = orientation;
        pixelListsMatrix{fiberID, frame} = pixelList;
        majorAxisLengthsMatrix(fiberID, frame) = majorAxisLength;
        minorAxisLengthsMatrix(fiberID, frame) = minorAxisLength;
    end
end

%% correct angles
cumulativeOrientationsMatrix = orientationsMatrix;
numFibers = size(orientationsMatrix, 1);
numFrames = size(orientationsMatrix, 2);

% Initialize cumulative angles with the first frame
for i = 1:numFibers
    for frame = 2:numFrames
        % Skip if current or previous frame is NaN
        if isnan(orientationsMatrix(i, frame)) || isnan(orientationsMatrix(i, frame - 1))
            continue;
        end
        % Calculate the angle difference
        angleDifference = orientationsMatrix(i, frame) - orientationsMatrix(i, frame - 1);
        
        % Correct for flipping
        if abs(angleDifference) > 90
            if angleDifference > 0
                angleDifference = angleDifference - 180;
            else
                angleDifference = angleDifference + 180;
            end
        end
        
        % Update the cumulative angle
        cumulativeOrientationsMatrix(i, frame) = cumulativeOrientationsMatrix(i, frame - 1) + angleDifference;
    end
end

orientationsMatrix = cumulativeOrientationsMatrix;

%% compute the end-points of the fiber
% for id = 1:size(positionsX,1)
%     for frame = 1:size(positionsX,2)
%         if ~isnan(positionsX(id,frame))
%             PixelList = pixelListsMatrix{id,frame};
%             % Calculate the pairwise distances between all points in PixelList
%             distances = pdist(PixelList);
%     
%             % Convert the distance vector to a square matrix
%             distMatrix = squareform(distances);
%     
%             % Find the pair of points with the maximum distance in the upper triangular part of the distance matrix
%             % Set the diagonal and lower triangular part to -inf to avoid selecting them
%             distMatrix = triu(distMatrix, 1);
%             [maxDist, maxIdx] = max(distMatrix(:));
%             
%             % Convert the linear index to row and column indices
%             [row, col] = ind2sub(size(distMatrix), maxIdx);
%             
%             % End-points of the fiber
%             EndPoints{id,frame}(1,:) =  PixelList(row, :);
%             EndPoints{id,frame}(2,:) =  PixelList(col, :);
%     
%             % % check
%     %         figure(1);clf;
%     %         plot(PixelList(:,2),PixelList(:,1),'k.','MarkerSize',10); daspect([1 1 1]); set(gcf,'Renderer','painters'); hold on;
%     %         plot(EndPoints{id,frame}(1,2),EndPoints{id,frame}(1,1),'r.','MarkerSize',30)
%     %         plot(EndPoints{id,frame}(2,2),EndPoints{id,frame}(2,1),'r.','MarkerSize',30)
%         end
% 
%     end
% end


%% save the data
% Save the matrices to files
save([p.main_folder,'only_fibers/','tracked_fibers.mat'], ...
    'positionsX','positionsY','orientationsMatrix',...
    'pixelListsMatrix',...
    'majorAxisLengthsMatrix','minorAxisLengthsMatrix','p');


%% Display results with random colors for each trajectory
if p.plot_images==1
    figure; set(gcf,'Position',[1 1 1200 900])
    hold on;
    colormap("gray")
    rng('default'); % For reproducibility, remove this line if you want different colors each run
    colors = rand(max(fiberIDs), 3); % Generate random colors

    for frame = 1:numFrames
        imagesc(readimage(imds,frame)); daspect([1 1 1]);
        hold on;
        for i = 1:size(trackedCentroids{frame}, 1)
            fiberID = trackedCentroids{frame}(i, 1);
            centroid = trackedCentroids{frame}(i, 2:3);
            orientation = trackedOrientations{frame}(i, 2);
            
            % Plot the trajectory
            if size(trajectories{fiberID}, 1) > 1
                plot(trajectories{fiberID}(:,1), trajectories{fiberID}(:,2), 'Color', colors(fiberID, :));
            end
            
            plot(centroid(1), centroid(2), 'o', 'Color', colors(fiberID, :), 'MarkerFaceColor', colors(fiberID, :));
            
        end
        hold off;
        title(['Frame ' num2str(frame)]);
        drawnow
    end
    hold off;
end

end