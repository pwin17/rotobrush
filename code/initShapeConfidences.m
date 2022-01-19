function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, MaskOutline, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.

    ShapeConfidences(1) = struct;
    s = size(LocalWindows);
    ShapeConfidences.Maps(s(1)) = struct;


    for i = 1:s(1)
        if (0<= ColorConfidences{i}) && (ColorConfidences{i}<=fcutoff)
            sigma = SigmaMin;
        else
            fcfcutoff = (fc-fcutoff)^R;
            sigma = SigmaMin + (A*fcfcutoff);
        end

        map = zeros(WindowWidth, WindowWidth);
        window_center = LocalWindows(i,:);
        window_center_x = window_center(1);
        window_center_y = window_center(2);

        offset = WindowWidth/2;
        x_offset0 = max(0,window_center_x-offset);
        y_offset0 = max(0,window_center_y-offset);
        
        for x = 1:WindowWidth
            for y = 1:WindowWidth
                x_r = x+x_offset0;
                y_r = y+y_offset0;
                d = distanceToMask(x_r,y_r,MaskOutline);

                map(y,x) = 1 - exp(-(d^2)/(sigma^2));
            end
        end

        ShapeConfidences.Maps(i).Map = map;

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