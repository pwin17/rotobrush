function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% TODO
%     opticFlow = opticalFlowHS;
    opticFlow = opticalFlowFarneback;
    warped_gray = rgb2gray(WarpedPrevFrame);
    current_gray = rgb2gray(CurrentFrame);
    flow = estimateFlow(opticFlow,warped_gray);
    flow = estimateFlow(opticFlow,current_gray);
    
    vx_avg = 0;
    vy_avg = 0;
    count = 0;
    
    for y = 1:height(Mask)
        for x = 1:width(Mask)
            if Mask(y,x)==true
                vx_avg = vx_avg + flow.Vx(y,x);
                vy_avg = vy_avg + flow.Vy(y,x);
                count = count + 1;
            end
        end
    end
    
    vx_avg = vx_avg/count;
    vy_avg = vy_avg/count;
    while abs(vx_avg) < 0.1
        vx_avg = vx_avg*10;
    end
    while abs(vy_avg) < 0.1
        vy_avg = vy_avg *10;
    end
%     h = figure;
%     movegui(h);
%     hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
%     hPlot = axes(hViewPanel);
%     plot(flow,'DecimationFactor',[5 5],'ScaleFactor',60,'Parent',hPlot);
    
%     vx_mean = mean(flow.Vx,'all');
%     vy_mean = mean(flow.Vy,'all');
    NewLocalWindows = [];
    for i = 1: length(LocalWindows)
       new_x = ceil(LocalWindows(i,1)+vx_avg*0.9);
       new_y = ceil(LocalWindows(i,2)+vy_avg);
       NewLocalWindows = [NewLocalWindows; new_x new_y];
    end
end

