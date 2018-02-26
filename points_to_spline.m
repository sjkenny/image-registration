% generate smoothing spline interpolation from hull points
% points = [x y] Nx2 matrix
% hull_idx is connectivity across points
% smoothingparam = [0,1] - lower is smoother
% output is x,y values of spline fit

function [xxx,yyy] = points_to_spline(hull_idx,points,smoothingparam)

extra = 0.2;   %add on to end of spline for continuity


hull_points_1 = points(hull_idx,:);
%add last point onto start of list for closed curve
hull_points_11 = cat(1,hull_points_1(end,:),hull_points_1);

t1 = 1:length(hull_points_11);

f = fit(t1',hull_points_11(:,1),'smoothingspline','SmoothingParam',smoothingparam);
f2 = fit(t1',hull_points_11(:,2),'smoothingspline','SmoothingParam',smoothingparam);
xxx = feval(f,linspace(1,length(t1)+extra,1000));
yyy = feval(f2,linspace(1,length(t1)+extra,1000));


plot(points(:,1),points(:,2),'k+')
hold on
plot(xxx,yyy,'c.','LineWidth',2)