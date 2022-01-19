function [mask, LocalWindows, newColorModels, ShapeConfidences] = ...
updateModels(...
    NewLocalWindows, ... %Centers of windows (x,y) based after optical flow
    LocalWindows, ... %LocalWindows holds x,y of all wanted local window centers
    CurrentFrame, ... %Current colored image frame
    warpedMask, ... %Logical warped mask
    warpedMaskOutline, ... %Logical warped mask outline
    WindowWidth, ... %width of window
    ColorModels, ... %Cell array of all local windows, each with a respective color confidence (ColorModels.Confidence)(Just a number)
    ShapeConfidences, ... %Struct with ShapeConfidences.Maps. That holds all wanted 40x40 local windows. 40x40 is logical
    ProbMaskThreshold, ... % >0.75 == foreground. >0.25==background
    fcutoff, ... %random value
    SigmaMin, ... %random value
    R, ... %2
    A ... %(SigmaMax - SigmaMin)/(1-fcutoff)^R
)
% UPDATEMODELS: update shape and color models, and apply the result to generate a new mask.
% Feel free to redefine this as several different functions if you prefer.

sc_size = size(ShapeConfidences.Maps);
newColorModels = struct;
newColorModels.Confidences = {}; % convert back later struct2cell
newColorModels.GMM = {};
newColorModels.pcx(sc_size(2)) = struct;
% not sure if I need to create a new shape Model
newShapeConfidences = ShapeConfidences;
disp('checking and calculating new shape and color models');
for i = 1: sc_size(2)
%     disp(i);
    map = newShapeConfidences.Maps(i).Map;
    map_size = size(map);
    window_center = NewLocalWindows(i,:);
    window_center_x = window_center(1);
    window_center_y = window_center(2);
    offset = WindowWidth/2;
    x_offset0 = max(0,window_center_x-offset);
    x_offset1 = min(width(CurrentFrame),window_center_x+offset);
    y_offset0 = max(0,window_center_y-offset);
    y_offset1 = min(height(CurrentFrame),window_center_y+offset);
    fg = [];
    bg = [];
    for j = 1:map_size(1)
       for k = 1:map_size(2)
          if map(j,k) > 0.75
              color = rgb2lab([CurrentFrame(j+y_offset0,k+x_offset0,1) CurrentFrame(j+y_offset0,k+x_offset0,2) CurrentFrame(j+y_offset0,k+x_offset0,3)]);
              fg = [fg; color];
          elseif map(j,k) > 0.25
              color = rgb2lab([CurrentFrame(j+y_offset0,k+x_offset0,1) CurrentFrame(j+y_offset0,k+x_offset0,2) CurrentFrame(j+y_offset0,k+x_offset0,3)]);
              bg = [bg ; color];
          end
       end
    end
    if length(fg)< 3
        inside_GMM = ColorModels.GMM{1,i}.foreground;
    else
        try
            inside_GMM = fitgmdist(fg, 3);
        catch exception
            rng(3); % Reset seed for common start values
            inside_GMM = fitgmdist(fg,3,'RegularizationValue',0.2);
        end
    end
    rng(1);
    if length(bg) < 3
        outside_GMM = ColorModels.GMM{1,i}.background;
    else
        try
            outside_GMM = fitgmdist(bg, 3);
        catch exception
            rng(3); % Reset seed for common start values
            outside_GMM = fitgmdist(bg,3,'RegularizationValue',0.2);
        end
    end
    total_pcx_new = 0;
    total_pcx_old = 0;
    oldGMM = ColorModels.GMM{1,i};
    old_fgGMM = oldGMM.foreground;
    old_bgGMM = oldGMM.background;
    for y = y_offset0:y_offset1
        for x = x_offset0:x_offset1
            color = rgb2lab([CurrentFrame(y,x,1) CurrentFrame(y,x,2) CurrentFrame(y,x,3)]);
            foreground = 0;
            if warpedMask(y,x) == 1
                foreground = 1;
            end
            pcf_new = posterior(inside_GMM, color);
            pcg_new = posterior(outside_GMM, color);
            pcx_new = pcf_new / (pcf_new + pcg_new);
            total_pcx_new = total_pcx_new + pcx_new;


            pcf_old = posterior(old_fgGMM, color);
            pcg_old = posterior(old_bgGMM, color);
            pcx_old = pcf_old / (pcf_old + pcg_old);
            total_pcx_old = total_pcx_old + pcx_old;

        end
    end

    map = zeros(WindowWidth, WindowWidth);
    if total_pcx_new >= total_pcx_old
%         disp('old color model');
        newColorModels.GMM{i} = struct;
        newColorModels.GMM{i}.foreground = old_fgGMM;
        newColorModels.GMM{i}.background = old_bgGMM;
        for y = y_offset0:y_offset1
            for x = x_offset0:x_offset1
                color = rgb2lab([CurrentFrame(y,x,1) CurrentFrame(y,x,2) CurrentFrame(y,x,3)]);
                foreground = 0;
                if warpedMask(y,x) == 1
                    foreground = 1;
                end
                pcf_new = posterior(old_fgGMM, color);
                pcg_new = posterior(old_bgGMM, color);
                pcx_new = pcf_new / (pcf_new + pcg_new);
                map(y-y_offset0+1,x-x_offset0+1) = pcx_new;
                d = distanceToMask(x,y,warpedMaskOutline);
                omega = exp(-(d^2)/(WindowWidth/2)^2);
                fc_numerator = abs(foreground - pcx_new) * omega;
                fc_denominator = omega;
            end
        end
        fc = 1 - (fc_numerator/fc_denominator);

    else
%         disp('new color model');
        newColorModels.GMM{i}.foreground = inside_GMM;
        newColorModels.GMM{i}.background = outside_GMM;
        for y = y_offset0:y_offset1
            for x = x_offset0:x_offset1
                color = rgb2lab([CurrentFrame(y,x,1) CurrentFrame(y,x,2) CurrentFrame(y,x,3)]);
                foreground = 0;
                if warpedMask(y,x) == 1
                    foreground = 1;
                end
                pcf_new = posterior(inside_GMM, color);
                pcg_new = posterior(outside_GMM, color);
                pcx_new = pcf_new / (pcf_new + pcg_new);
                map(y-y_offset0+1,x-x_offset0+1) = pcx_new;
                d = distanceToMask(x,y,warpedMaskOutline);
                omega = exp(-(d^2)/(WindowWidth/2)^2);
                fc_numerator = abs(foreground - pcx_new) * omega;
                fc_denominator = omega;
            end
        end
        fc = 1 - (fc_numerator/fc_denominator);
        if (0<= fc) && (fc<=fcutoff)
            sigma = SigmaMin;
        else
            fcfcutoff = (fc-fcutoff)^R;
            sigma = SigmaMin + (A*fcfcutoff);
        end
        shape_map = zeros(WindowWidth, WindowWidth);
        window_center = NewLocalWindows(i,:);
        window_center_x = window_center(1);
        window_center_y = window_center(2);

        offset = WindowWidth/2;
        x_offset0 = max(0,window_center_x-offset);
        y_offset0 = max(0,window_center_y-offset);

        for x = 1:WindowWidth
            for y = 1:WindowWidth
                x_r = x+x_offset0;
                y_r = y+y_offset0;
                d = distanceToMask(x_r,y_r,warpedMaskOutline);

                shape_map(y,x) = 1 - exp(-(d^2)/(sigma^2));
            end
        end

        newShapeConfidences.Maps(i).Map = shape_map;
    end
    newColorModels.Confidences{i} = fc;
    newColorModels.pcx(i).Map = map;
end

s_nlw = size(NewLocalWindows);
newColorModels.ForegroundProb(s_nlw(1)) = struct;

disp('Calculating Foreground Probs');
for i = 1:s_nlw(1)
    map = zeros(WindowWidth, WindowWidth);
    window_center = NewLocalWindows(i,:);
    window_center_x = window_center(1);
    window_center_y = window_center(2);

    offset = WindowWidth/2;
    x_offset0 = max(0,window_center_x-offset);
    y_offset0 = max(0,window_center_y-offset);
    shapeMap = newShapeConfidences.Maps(i).Map;
    pcx = newColorModels.pcx(i).Map;
    for y = 1:WindowWidth
        for x = 1:WindowWidth
           x_r = x+x_offset0;
           y_r = y+y_offset0;
           if warpedMask(y_r,x_r) == 1
               fg = 1;
           else
               fg = 0;
           end
           pkfx = (shapeMap(y,x)*fg) + ((1-shapeMap(y,x))*pcx(y,x));
           map(y,x) = pkfx;
        end
    end
    newColorModels.ForegroundProb(i).Map = map; %pFkx 40x40
end

ShapeConfidences = newShapeConfidences;
ColorModels = newColorModels;

RealValuedProbabilityMap(1) = struct;
w = width(CurrentFrame);
h = height(CurrentFrame); % change this to take height of image since different input sets have different dimensions

RealValuedProbabilityMap.Map= zeros(h, w); %Creates proper amount of new local windows
NewForegroundProbability = 0;

%%%%Imagining p_fx is the output of part 6
disp('Merging Local Windows');
%Go through every pixel in frame
for y = 1:h
    for x = 1:w
        %Go through all local windows to check if current pixel is in a
        %window
        NewForegroundProbability = 0; %resets value per pixel
        numerator = 0;
        denominator = 0;
        isInWindow = false;
        for k = 1:length(NewLocalWindows)
            if x > NewLocalWindows(k,1) - round(WindowWidth/2) && x < NewLocalWindows(k,1) + round(WindowWidth/2)
               if y > NewLocalWindows(k,2) - round(WindowWidth/2) && y < NewLocalWindows(k,2) + round(WindowWidth/2)
                    isInWindow = true;
                    topLeftCornerStart_x = NewLocalWindows(k,1) - round(WindowWidth/2); %x indice of topleft corner of current local window
                    topLeftCornerStart_y = NewLocalWindows(k,2) - round(WindowWidth/2); %y indice of topleft corner of current local window
                    %Move topLeftCorner_x/y accordingly to land on current
                    %pixel x,y and find the right value in the current
                    %local window k.
                    p_fXCoord = 1; %matrix index x to get p_f(x) of current local window
                    p_fYCoord = 1; %matrix index y to get p_f(x) of current local window
                    %Get matrix index x for current local window
                    while topLeftCornerStart_x < x
                       topLeftCornerStart_x = topLeftCornerStart_x + 1; 
                       p_fXCoord = p_fXCoord + 1;
                    end
                    %Get matrix index y for current local window
                    while topLeftCornerStart_y < y
                       topLeftCornerStart_y = topLeftCornerStart_y + 1; 
                       p_fYCoord = p_fYCoord + 1;
                    end
                    %%Is the indexing for p_fx.Maps(k) right??
                    %%ColorModels.ForegroundProb.Map(k).Map(
%                     disp(ColorModels.ForegroundProb(k).Map(p_fXCoord, p_fYCoord));
%                     disp(p_fXCoord);
%                     disp(p_fYCoord);
                    numerator = numerator + ColorModels.ForegroundProb(k).Map(p_fYCoord,p_fXCoord) * (getDist(x,y,p_fXCoord,p_fYCoord) + 0.1)^(-1);
                    denominator = denominator + (getDist(x,y,p_fXCoord,p_fYCoord) + 0.1)^(-1);
               end
            end
        end

        %If the pixel isn't in a local window
        % it should be the original mask value because the windows are only
        % in the boundary
        if isInWindow == false
            RealValuedProbabilityMap.Map(y,x) = warpedMask(y,x); 
        else %if the pixel is in a local window
            NewForegroundProbability = numerator/denominator;
            RealValuedProbabilityMap.Map(y,x) = NewForegroundProbability;
        end
    end


end

mask = RealValuedProbabilityMap.Map;
% mask_after = RealValuedProbabilityMap.Map;
LocalWindows = NewLocalWindows; %check

disp('Creating new Mask')

% L = superpixels(CurrentFrame,4000);
%Loop through mask now. If the foreground probability for a pixel is
%greater than established threshold (ProbMaskThreshold), then that pixel
%will be part of the foreground. Change that probability to 1. If below
%threshold, it's background so change probability to 0.
for y = 1:h
    for x = 1:w
        if mask(y,x) >= ProbMaskThreshold
            mask(y,x) = 1;
        else
            mask(y,x) = 0;
        end
    end
end

% background_mask = ~(mask);
% mask = lazysnapping(CurrentFrame,L,mask,background_mask);


end
function distance = getDist(x1, y1, x2, y2)
    distance = sqrt((x2-x1)^2 + (y2-y1)^2);
end
function dist = distanceToMask(X,Y,mask)
    s = size(mask);
    dist = Inf;
    for y = 1:s(1)
        for x = 1:s(2)
            if mask(y,x) == 1
                d = (Y - y)^2 + (X - x)^2;
                d = sqrt(d);
                if d < dist
                    dist = d;
                end
            end
        end
    end
end