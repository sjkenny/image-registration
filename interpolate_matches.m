% interpolate and quantify all Brp coords
% load cleaned matches (match_idx_crop from quantify_matches)
% load matches struct
% for cleaned match output from reg_clustering_affine_v3, use match_idx
%%


storm_centers = out.storm_centers;
fcn_centers = out.fcn_centers;
match_idx = out.match_idx;

fcn_match = fcn_centers(match_idx(:,2),:);
storm_match = storm_centers(match_idx(:,1),:);

shift = fcn_match-storm_match;

xFit = fit([storm_match],shift(:,1),'thinplateinterp');
yFit = fit([storm_match],shift(:,2),'thinplateinterp');

x_add = feval(xFit,storm_centers);
y_add = feval(yFit,storm_centers);

fcn_center_predict = bsxfun(@plus,storm_centers,[x_add y_add]);
%%
% plot(storm_centers(:,1),storm_centers(:,2),'k.')
% hold on
% axis equal
% plot(fcn_center_predict(:,1),fcn_center_predict(:,2),'m.')
% 
% for j=1:length(storm_centers)
% 
%     x=[fcn_center_predict(j,1) storm_centers(j,1)];
%     y=[fcn_center_predict(j,2) storm_centers(j,2)];
% 
%     line(x,y)
% 
% end


save('points_interp','fcn_center_predict','storm_centers')

[x,y] = meshgrid(linspace(min(storm_centers(:,1)),max(storm_centers(:,1)),100),...
            linspace(min(storm_centers(:,2)),max(storm_centers(:,2)),100));

u = feval(xFit,x,y);
v = feval(yFit,x,y);
% quiver(x,y,u,v,'Color','g','AutoScale','off')
hold on
plot(storm_centers(:,1),storm_centers(:,2),'k.')
hold on
axis equal
plot(fcn_center_predict(:,1),fcn_center_predict(:,2),'m.')
% interpolated matches

for j=1:length(storm_centers)

    x1=[fcn_center_predict(j,1) storm_centers(j,1)];
    y1=[fcn_center_predict(j,2) storm_centers(j,2)];

    line(x1,y1)
    

end
% matches
for j=1:length(match_idx)
    x=[fcn_centers(match_idx(j,2),1) storm_centers(match_idx(j,1),1)];
    y=[fcn_centers(match_idx(j,2),2) storm_centers(match_idx(j,1),2)];

    line(x,y,'color','r')
end
    