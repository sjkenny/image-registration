%% storm overlay

addpath ../common 

if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.png',LastFolder);
else
    GetFileName='*.png';
end


[FileNameL,PathNameL] = uigetfile(GetFileName,'Select Image 1');

GetFileName=sprintf('%s/*.ini',PathNameL);
[FileNameInf,PathNameInf] = uigetfile(GetFileName,'Select inf file');

ini_name = sprintf('%s%s',PathNameInf,FileNameInf);


a = ReadInfoFile(ini_name);
zoom = a.Zoom;
xStart = a.DispX-0.5;
yStart = a.DispY-0.5;

LeftFile = sprintf('%s%s',PathNameL,FileNameL);
LeftImg=imread(LeftFile);
%% plot hull over storm img
fcn_use = [reg_evoked_Ib; reg_evoked_Is];
fcn_use_scale = bsxfun(@minus,fcn_use,[xStart yStart]);
fcn_use_scale(:,1) = fcn_use_scale(:,1)+150;
fcn_use_scale = fcn_use_scale.*zoom;
% brp_use_scale = brp_use;
% brp_use_scale_1 = brp_use1;

% brp_use_scale = bsxfun(@minus,brp_use_scale,[xStart yStart]);
% % brp_use_scale_1 = bsxfun(@minus,brp_use_scale_1,[xStart yStart]);
% 
% brp_use_scale = brp_use_scale*zoom;
% brp_use_scale_1 = brp_use_scale_1*zoom;

LeftImgPad = padarray(LeftImg,[0 450],'post');
clf
imshow(LeftImgPad)
hold on
% plot_hull(brp_use_scale,hull_idx);
% plot_hull(brp_use_scale_1,hull_idx1);
plot(fcn_use_scale(:,1),fcn_use_scale(:,2),'c.','MarkerSize',0.5)
%%
imshow(LeftImg)
hold on
% plot_hull(brp_use_scale,hull_idx);
% plot_hull(brp_use_scale_1,hull_idx1);
plot(xxx,yyy,'c','LineWidth',2)

%% test spline fitting
% get hull points and indices from compute_mask

extra = 0.2;   %add on to end of spline for continuity

smoothingparam = 0.99;

hull_points_1 = brp_use_scale(hull_idx,:);
%add last point onto start of list for closed curve
hull_points_11 = cat(1,hull_points_1(end,:),hull_points_1);

t1 = 1:length(hull_points_11);

f = fit(t1',hull_points_11(:,1),'smoothingspline','SmoothingParam',smoothingparam);
f2 = fit(t1',hull_points_11(:,2),'smoothingspline','SmoothingParam',smoothingparam);
xxx = feval(f,linspace(1,length(t1)+extra,100));
yyy = feval(f2,linspace(1,length(t1)+extra,100));
% 
% hull_points_2 = brp_use_scale_1(hull_idx1,:);
%add last point onto start of list for closed curve
% hull_points_22 = cat(1,hull_points_2(end,:),hull_points_2);

% t1 = 1:length(hull_points_22);

% f = fit(t1',hull_points_22(:,1),'smoothingspline','SmoothingParam',smoothingparam);
% f2 = fit(t1',hull_points_22(:,2),'smoothingspline','SmoothingParam',smoothingparam);
% x2 = feval(f,linspace(1,length(t1)+extra,100));
% y2 = feval(f2,linspace(1,length(t1)+extra,100));


imshow(LeftImg)
hold on
% plot_hull(brp_use_scale,hull_idx);
% plot_hull(brp_use_scale_1,hull_idx1);
plot(xxx,yyy,'c','LineWidth',2)
% plot(x2,y2,'y','LineWidth',2)

%% plot matches over storm img


match_idx = out.match_idx;

brp_use = out.storm_centers(match_idx(:,1),:);
fcn_use = out.fcn_centers(match_idx(:,2),:);

brp_all = storm_centers_brp_include;
fcn_all = out.fcn_centers;
fcn_all(:,1) = fcn_all(:,1)+256;
brp_all_scale = bsxfun(@minus,brp_all,shift);
fcn_all_scale = bsxfun(@minus,fcn_all,shift);

% fcn_center_predict_scale = fcn_center_predict;
fcn_center_scale = fcn_use;
fcn_center_scale(:,1) = fcn_center_scale(:,1)+256;
storm_centers_scale = brp_use;
% fcn_center_predict_scale(:,1) = fcn_center_predict_scale(:,1)+256;
%coords start halfway through pixel

shift = [xStart-0.5 yStart-0.5];


% storm_centers = out.storm_centers;
% fcn_centers = fcn_centers.*zoom;

storm_centers_scale = bsxfun(@minus,storm_centers_scale,shift);
% fcn_center_predict_scale = bsxfun(@minus,fcn_center_predict_scale,shift);
fcn_center_scale = bsxfun(@minus,fcn_center_scale,shift);
storm_centers_scale = storm_centers_scale.*zoom;
fcn_center_scale = fcn_center_scale.*zoom;

fcn_all_scale = fcn_all_scale.*zoom;
brp_all_scale = brp_all_scale.*zoom;
% fcn_center_predict_scale = fcn_center_predict_scale.*2;

imshow(LeftImg)

plot(brp_all_scale(:,1),brp_all_scale(:,2),'yo')
plot(fcn_all_scale(:,1),fcn_all_scale(:,2),'yo')
plot(storm_centers_scale(:,1),storm_centers_scale(:,2),'wo')
plot(fcn_center_scale(:,1),fcn_center_scale(:,2),'wo')

%%
for j=1:length(brp_use)
    
    x=[fcn_center_scale(j,1) storm_centers_scale(j,1)];
    y=[fcn_center_scale(j,2) storm_centers_scale(j,2)];

    line(x,y)

end


