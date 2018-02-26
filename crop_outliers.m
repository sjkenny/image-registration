function [idx_include] = crop_outliers(points)


KeepGoing=1;
idx_include = 1:length(points);
outlier_idx = [];
while KeepGoing==1
    
    clf
    plot(points(idx_include,1),points(idx_include,2),'b+','MarkerSize',5)
    axis equal
    ax=gca;
    try
        rect=getrect(ax);
    catch
        KeepGoing=0;
        break
    end
    % find points within selected rectangle
    outlier_idx_now = find(points(:,1)>rect(1)&points(:,2)>rect(2)&...
    points(:,1)<(rect(1)+rect(3))&points(:,2)<(rect(2)+rect(4)));
    idx_include = 1:length(points);
    outlier_idx = cat(1,outlier_idx,outlier_idx_now);   %add to list
%     idx_storm_include = 1:length(brp_all); 
    idx_include(intersect(idx_include,outlier_idx))=[];
    clf
    plot(points(idx_include,1),points(idx_include,2),'b+','MarkerSize',5)
    axis equal
end

