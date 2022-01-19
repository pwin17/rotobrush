% MyRotobrush.m  - UMD CMSC426, Fall 2018
% This is the main script of your rotobrush project.
% We've included an outline of what you should be doing, and some helful visualizations.
% However many of the most important functions are left for you to implement.
% Feel free to modify this code as you see fit.

% Some parameters you need to tune:
WindowWidth = 40; % how big the feature descriptors are?  
ProbMaskThreshold = 0.65; % probability threshold for being part of the image?
NumWindows= 25; % number of equally spaced feature descriptors? 27 for frames 3
BoundaryWidth = 3; % how much do we want to separate between the inner and outer region

% Load images:
fpath = '/Users/pyone/Desktop/Fall2021/vision/CMSC426_P3/input'; % CHANGE THIS
outpath = '/Users/pyone/Desktop/Fall2021/vision/CMSC426_P3/output/f1/';
files = dir(fullfile(fpath, '*.jpg'));
imageNames = zeros(length(files)-1,1);
images = cell(length(files),1);

for i=1:length(files)-1
    imageNames(i) = str2double(strtok(files(i).name,'.jpg'));
end

imageNames = sort(imageNames);
imageNames = num2str(imageNames);
imageNames = strcat(imageNames, '.jpg');

for i=1:length(files)-1
    images{i} = im2double(imread(fullfile(fpath, strip(imageNames(i,:)))));
end

mask = imread(strcat(fpath,'/Mask1.jpg'));
mask = imbinarize(mask);
gt_mask = mask;
checkImg = imoverlay(images{1}, boundarymask(mask,8), 'red');
imshow(checkImg);
imwrite(checkImg,strcat(outpath,'f1_overlay.jpg'));
set(gca,'position',[0 0 1 1],'units','normalized');
F = getframe(gcf);
[I,~] = frame2im(F);

imwrite(I, fullfile(outpath, strip(imageNames(1,:))));
outputVideo = VideoWriter(fullfile(outpath,'video1.mp4'),'MPEG-4');
open(outputVideo);
writeVideo(outputVideo,I);


% Sample local windows and initialize shape+color models:
[mask_outline, LocalWindows] = initLocalWindows(images{1},mask,NumWindows,WindowWidth,true);

% increase the boundary width
img_size = size(images{1});
mask_outline = increase_boundary(mask_outline,BoundaryWidth,img_size);

ColorModels = ...
    initColorModels(images{1},mask,mask_outline,LocalWindows,BoundaryWidth,WindowWidth);

% You should set these parameters yourself:
fcutoff = 0.85; % value from paper
SigmaMin = 2; % value from paper
SigmaMax = WindowWidth; % max from paper
R = 2; % parameter to calculate sigma
Anumer = SigmaMax-SigmaMin;
Adenom = (1-fcutoff)^R;
A = Anumer/Adenom; % parameter to calculate sigma
ShapeConfidences = ...
    initShapeConfidences(LocalWindows,ColorModels.Confidences, mask_outline,...
    WindowWidth, SigmaMin, A, fcutoff, R);

% Show initial local windows and output of the color model:
imshow(images{1})
hold on
showLocalWindows(LocalWindows,WindowWidth,'r.');
hold off

set(gca,'position',[0 0 1 1],'units','normalized')
F = getframe(gcf);
[I,~] = frame2im(F);

hold on
showColorConfidences(images{1},mask_outline,ColorModels.Confidences,LocalWindows,WindowWidth);
hold off
% saveas(gca,strcat(outpath,'f1_colorconfidences.jpg'));
%%% MAIN LOOP %%%
% Process each frame in the video.
for prev=1:(length(files)-2)
    curr = prev+1;
    fprintf('Current frame: %i\n', curr)
    try
        %%% Global affine transform between previous and current frames:
        [warpedFrame, warpedMask, warpedMaskOutline, warpedLocalWindows] = calculateGlobalAffine(images{prev}, images{curr}, mask, LocalWindows, BoundaryWidth);

        %%% Calculate and apply local warping based on optical flow:
        NewLocalWindows = ...
            localFlowWarp(warpedFrame, images{curr}, warpedLocalWindows,warpedMask,WindowWidth);

        % Show windows before and after optical flow-based warp:
        imshow(images{curr});
        hold on
        showLocalWindows(warpedLocalWindows,WindowWidth,'r.');
        showLocalWindows(NewLocalWindows,WindowWidth,'b.');
        hold off

        %%% UPDATE SHAPE AND COLOR MODELS:
        % This is where most things happen.
        % Feel free to redefine this as several different functions if you prefer.
        [ ...
            mask, ...
            LocalWindows, ...
            ColorModels, ...
            ShapeConfidences, ...
        ] = ...
        updateModels(...
            NewLocalWindows, ...
            LocalWindows, ...
            images{curr}, ...
            warpedMask, ...
            warpedMaskOutline, ...
            WindowWidth, ...
            ColorModels, ...
            ShapeConfidences, ...
            ProbMaskThreshold, ...
            fcutoff, ...
            SigmaMin, ...
            R, ...
            A ...
        );

        mask_outline = bwperim(mask,4);
        mask_outline = increase_boundary(mask_outline,BoundaryWidth,size(images{curr}));
        
        % Write video frame:
        imshow(imoverlay(images{curr}, boundarymask(mask,8), 'red'));
        imwrite(imoverlay(images{curr}, boundarymask(mask,8), 'red'),strcat(outpath,'f',num2str(curr),'_overlay.jpg'));
        set(gca,'position',[0 0 1 1],'units','normalized')
        F = getframe(gcf);
        [I,~] = frame2im(F);
        imwrite(I, fullfile(outpath, strip(imageNames(curr,:))));
        writeVideo(outputVideo,I);

        imshow(images{curr})
        hold on
        showLocalWindows(LocalWindows,WindowWidth,'r.');
        hold off
        set(gca,'position',[0 0 1 1],'units','normalized')
        F = getframe(gcf);
        [I,~] = frame2im(F);
    catch exception
        imshow(imoverlay(images{curr}, boundarymask(mask,8), 'red'));
        imwrite(imoverlay(images{curr}, boundarymask(mask,8), 'red'),strcat(outpath,'f',num2str(curr),'_overlay.jpg'));
        set(gca,'position',[0 0 1 1],'units','normalized')
        F = getframe(gcf);
        [I,~] = frame2im(F);
        imwrite(I, fullfile(outpath, strip(imageNames(curr,:))));
        writeVideo(outputVideo,I);
    end
end

close(outputVideo);

function new_mask_outline = increase_boundary(mask_outline,boundary_width,img_size)
    boundary_offset = fix(boundary_width/2);
    indexes = [];
    for y = 1:img_size(1)
        for x = 1:img_size(2)
            if mask_outline(y,x) == 1
                indexes = [indexes; y x];
            end
        end
    end

    for i = 1: length(indexes)
        yx = indexes(i,:);
        for b_i = 1: boundary_offset
            mask_outline(yx(1)-b_i,yx(2)) = 1;
            mask_outline(yx(1)+b_i,yx(2)) = 1;
            mask_outline(yx(1),yx(2)-b_i) = 1;
            mask_outline(yx(1),yx(2)+b_i) = 1;
        end
    end
    new_mask_outline = mask_outline;
end