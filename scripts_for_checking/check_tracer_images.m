clc; clear; close all
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath('C:\Users\Corsair\Desktop\Vlad\___common\')


%% input
main_folder = 'E:\ML_paper\Re1000_fiber_v4\Re1000_Fiber_3mm_01\loop=0\';


%% load images
% image data store
imds = imageDatastore([main_folder,'only_tracers/'],'FileExtensions','.tif');

for i = 1:size(imds.Files,1)
    
    figure(1);clf; box on; daspect([1 1 1]); set(gcf,'WindowState', 'maximized');
    imagesc(readimage(imds,i));
    daspect([1 1 1]);
    clim([0 4095])
    colorbar
    title('Tracers');
    pause(0.1)

end