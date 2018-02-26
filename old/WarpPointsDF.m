% Warp coordinates using displacement field
% pixel center should be +.5 in each dimension
% doesn't work for negative coordinates

function [out_x,out_y] = WarpPointsDF(x,y,dField)

out_x = x.*0;
out_y = y.*0;

DX = dField(:,:,1);
DY = dField(:,:,2);


xShift = x-0.5;
yShift = y-0.5;

xHigh = ceil(xShift);
xLow = floor(xShift);
yHigh = ceil(yShift);
yLow = floor(yShift);

for i=1:length(x)
    xNow = x(i);
    yNow = y(i);
    xLow_now = max(xLow(i),1);
    yLow_now = max(yLow(i),1);
                    % value*slope+intercept
    val_x_low = (xShift(i)-xLow_now)*(DX(yLow_now,xHigh(i))-DX(yLow_now,xLow_now))+DX(yLow_now,xLow_now);
    val_x_high = (xShift(i)-xLow_now)*(DX(yHigh(i),xHigh(i))-DX(yHigh(i),xLow_now))+DX(yHigh(i),xLow_now);
    %value of x shift
    val_x = (yShift(i)-yLow_now)*(val_x_high-val_x_low)+val_x_low;
    
    val_y_low = (yShift(i)-yLow_now)*(DY(yHigh(i),xLow_now)-DY(yLow_now,xLow_now))+DY(yLow_now,xLow_now);
    val_y_high = (yShift(i)-yLow_now)*(DY(yHigh(i),xHigh(i))-DY(yLow_now,xHigh(i)))+DY(yLow_now,xHigh(i));
    val_y = (xShift(i)-xLow_now)*(val_y_high-val_y_low)+val_y_low;
    
    out_x(i) = x(i)+val_x;
    out_y(i) = y(i)+val_y;
end
% plot(x,y,'k.')
% hold on
% plot(out_x,out_y,'m.')
    
    



