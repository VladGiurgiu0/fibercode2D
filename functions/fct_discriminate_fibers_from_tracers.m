function [] = fct_discriminate_fibers_from_tracers(p)
%% load images
% image data store
imds = imageDatastore([p.main_folder],'FileExtensions','.tif');

%% make folders
if ~exist([p.main_folder,'only_fibers/'],'dir'); mkdir([p.main_folder,'only_fibers/']); end
if ~exist([p.main_folder,'only_tracers/'],'dir'); mkdir([p.main_folder,'only_tracers/']); end

%% process images
fiberMask_cell = cell(size(imds.Files, 1), 1); % Proper initialization for use within parfor

parfor i = 1:size(imds.Files,1)
    Img = readimage(imds,i);
    
    % Normalize the image to the range [0, 1]
    normalizedImg = mat2gray(Img);
    
    % Thresholding to create a binary image
    binaryImg = imbinarize(normalizedImg,p.threshold_for_binarization_keep_fibers);
    
    % Measure properties of image regions
    stats = regionprops(binaryImg, 'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList','Centroid','Orientation');
    
    % Initialize masks for fibers
    fiberMask = false(size(Img));
    
    % Separate fibers from tracers based on aspect ratio
    for k = 1:length(stats)
        %aspectRatio = stats(k).MajorAxisLength / stats(k).MinorAxisLength;
        fiber_length = stats(k).MajorAxisLength;
        if fiber_length > p.length_threshold_keep_fibers % Adjust this threshold based on your criteria
            fiberMask(stats(k).PixelIdxList) = true;
        end
    end
    
    %fiber_Img = Img .* uint16(fiberMask);


    % Delete fibers from tracer images
    % Thresholding to create a binary image
    binaryImg_remove_fibers = imbinarize(normalizedImg,p.threshold_for_binarization_remove_fibers);
    
    % Measure properties of image regions
    stats_remove_fibers = regionprops(binaryImg_remove_fibers, 'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList','Centroid','Orientation');
    
    % Initialize masks for fibers
    fiberMask_remove_fibers = false(size(Img));
    
    % Separate fibers from tracers based on aspect ratio
    for k = 1:length(stats_remove_fibers)
        %aspectRatio = stats(k).MajorAxisLength / stats(k).MinorAxisLength;
        fiber_length = stats_remove_fibers(k).MajorAxisLength;
        if fiber_length > p.length_threshold_remove_fibers % Adjust this threshold based on your criteria
            fiberMask_remove_fibers(stats_remove_fibers(k).PixelIdxList) = true;
        end
    end

    se = strel('disk', p.dilation_size); % You can adjust the size and shape of the structuring element
    dilatedFiberMask_remove_fibers = imdilate(fiberMask_remove_fibers, se);
    tracer_Img = Img .* uint16(~dilatedFiberMask_remove_fibers);

    % Process fiber images to improve fiber image quality
    % deletes all fibers which touch the border
    %figure(1);clf;
    %subplot(1,2,1)
    %imagesc(fiber_Img); axis equal
    fiberMask = imclearborder(fiberMask); 
    %fiber_Img = fiber_Img.*uint16(fiberMask);
    %subplot(1,2,2)
    %imagesc(fiber_Img); axis equal

    % cleans the fiber edges
    %figure(1);clf;
    %subplot(1,3,1)
    %imagesc(fiber_Img); axis equal
    %subplot(1,3,2)
    fiberMask = bwmorph(fiberMask,'spur');
    %fiberMask = bwmorph(fiberMask,'majority');
    %imagesc(fiber_Img.*uint16(fiberMask)); axis equal
    
    %subplot(1,3,3)
    %fiber_Img = imgaussfilt(fiber_Img, 1); % Adjust the filter size as needed
    %imagesc(fiber_Img); axis equal

    fiber_Img =Img.*uint16(fiberMask);

    %fiber_Img = imdilate(fiber_Img, se);

    if p.plot_images==1
        % Display results
        %fig1=figure(1); fig1.Position=[1 1 1800 900]; clf;
        %subplot(3,1,1);
        figure(1); set(gcf,'Position',[100 1 1300 800]);
        imagesc(Img);
        daspect([1 1 1]);
        clim([0 4095])
        %colorbar
        xticklabels([])
        yticklabels([])
        title('Original Image');
        
        %subplot(3,1,2);
        figure(2); set(gcf,'Position',[200 1 1300 800]);
        imagesc(fiber_Img);
        daspect([1 1 1]);
        clim([0 4095])
        %colorbar
        xticklabels([])
        yticklabels([])
        title('Separated Fibers');
        
        %subplot(3,1,3);
        figure(3); set(gcf,'Position',[300 1 1300 800]);
        imagesc(tracer_Img);
        daspect([1 1 1]);
        clim([0 4095])
        %colorbar
        xticklabels([])
        yticklabels([])
        title('Separated Tracers');
        %pause(0.01)
        pause
    end

    % Save images
    imwrite(fiber_Img, [p.main_folder,'only_fibers/','fibers_',num2str(i,'%04d'),'.tif']);
    imwrite(tracer_Img, [p.main_folder,'only_tracers/','tracers_',num2str(i,'%04d'),'.tif']);

    fiberMask_cell{i} = fiberMask;

    disp(['Processed image number: ',num2str(i),' out of ', num2str(size(imds.Files,1))])

end

save([p.main_folder,'only_fibers/','mask_fibers.mat'],"fiberMask_cell",'-v7.3');

end