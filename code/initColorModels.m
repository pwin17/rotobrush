function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
% INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
%
% Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.

ColorModels(1) = struct;
s = size(LocalWindows);
ColorModels.Confidences = {}; % convert back later struct2cell
ColorModels.GMM = {};
ColorModels.pcx(s(1)) = struct;
for i = 1:s(1)
    disp(i);
    window_center = LocalWindows(i,:);
    window_center_x = window_center(1);
    window_center_y = window_center(2);
    
    offset = WindowWidth/2;
    x_offset0 = max(0,window_center_x-offset);
    x_offset1 = min(width(IMG),window_center_x+offset);
    y_offset0 = max(0,window_center_y-offset);
    y_offset1 = min(height(IMG),window_center_y+offset);
    
    colors_inside_mask = [];
    colors_outside_mask = [];

    for y = y_offset0:y_offset1
        for x = x_offset0:x_offset1
            color = rgb2lab([IMG(y,x,1) IMG(y,x,2) IMG(y,x,3)]);
            if Mask(y,x) == 1
                if MaskOutline(y,x) ~= 1
                    colors_inside_mask = [colors_inside_mask; color];
                end
            else
                if MaskOutline(y,x) ~= 1
                    colors_outside_mask = [colors_outside_mask; color];
                end
            end
        end
    end
    rng(1);
%     disp('inside GMM');
    try
        inside_GMM = fitgmdist(colors_inside_mask, 3);
    catch exception
        disp(exception);
        rng(3); % Reset seed for common start values
        inside_GMM = fitgmdist(colors_inside_mask,3,'RegularizationValue',0.2);
    end
%     disp('outside GMM');
    rng(1);
    try
        outside_GMM = fitgmdist(colors_outside_mask, 3);
    catch
        rng(3); % Reset seed for common start values
        outside_GMM = fitgmdist(colors_outside_mask,3,'RegularizationValue',0.2);
    end
    ColorModels.GMM{i} = struct;
    ColorModels.GMM{i}.foreground = inside_GMM;
    ColorModels.GMM{i}.background = outside_GMM;
    fc_numerator = 0;
    fc_denominator = 0;
    map = zeros(WindowWidth, WindowWidth);
    for y = y_offset0:y_offset1
        for x = x_offset0:x_offset1
            
            color = rgb2lab([IMG(y,x,1) IMG(y,x,2) IMG(y,x,3)]);
            foreground = 0;
            if Mask(y,x) == 1
                foreground = 1;
            end
            pcf = posterior(inside_GMM, color);
            pcg = posterior(outside_GMM, color);
            pcx = pcf / (pcf + pcg);
            if pcx > 0.5
                pcx = pcx+0.1;
            elseif pcx <0.5
                pcx = pcx-0.1;
            end
            pcx = 1-pcx;
            map(y-y_offset0+1,x-x_offset0+1) = pcx;
            d = distanceToMask(x,y,MaskOutline);
            omega = exp(-(d^2)/(WindowWidth/2)^2);
            fc_numerator = fc_numerator + (abs(foreground - pcx) * omega);
            fc_denominator = fc_denominator + omega;
        end
    end
    
    fc = 1 - (fc_numerator/fc_denominator);
    ColorModels.Confidences{i} = fc;
    ColorModels.pcx(i).Map = map;
    
end

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