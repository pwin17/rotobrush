function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows,boundary_width)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    im1_gray = rgb2gray(IMG1);
    im2_gray = rgb2gray(IMG2);
    
    im1_corners = detectHarrisFeatures(im1_gray);
    [im1_features,im1_valid_pts] = extractFeatures(im1_gray,im1_corners);
    im2_corners = detectHarrisFeatures(im2_gray);
    [im2_features,im2_valid_pts] = extractFeatures(im2_gray,im2_corners);
    
    index_pairs = matchFeatures(im1_features,im2_features);
    
    im1_matched_pts = im1_valid_pts(index_pairs(:,1));
    im2_matched_pts = im2_valid_pts(index_pairs(:,2));
    
    [tform,inlierIdx] = estimateGeometricTransform2D(im2_matched_pts,im1_matched_pts,'affine');

    outputView = imref2d(size(IMG1));
    WarpedFrame = imwarp(IMG1,tform,'OutputView',outputView);
    WarpedMask = imwarp(Mask,tform,'OutputView',outputView);
    WarpedMaskOutline = bwperim(WarpedMask,4);
    WarpedMaskOutline = increase_boundary(WarpedMaskOutline,boundary_width,size(WarpedMask));
    WarpedLocalWindows = [];
    for i = 1:length(Windows)
        new_coords = [Windows(i,1) Windows(i,2) 1] * tform.T;
        WarpedLocalWindows = [WarpedLocalWindows; new_coords(1) new_coords(2)];
    end
    
end
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
