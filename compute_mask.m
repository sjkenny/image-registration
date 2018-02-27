% test version
% compute inlier mask of triangulation
% SJK 1/18/18
%load in out struct

addpath ../common
brp_all = double(out.storm_centers);

% delete outliers
% keep track of included and deleted indices so matches can be referenced back to
% original list
% make this into separate function - crop_outliers

KeepGoing=1;
idx_storm_include = 1:length(brp_all);
storm_outlier_idx = [];
while KeepGoing==1
    
    clf
    plot(brp_all(idx_storm_include,1),brp_all(idx_storm_include,2),'b+','MarkerSize',5)
    hold on
    plot(out.fcn_loc_coords(:,1),out.fcn_loc_coords(:,2),'m.')
    axis equal
    ax=gca;
    try
        rect=getrect(ax);
    catch
        KeepGoing=0;
        break
    end
    
    outlier_idx_now = find(brp_all(:,1)>rect(1)&brp_all(:,2)>rect(2)&...
    brp_all(:,1)<(rect(1)+rect(3))&brp_all(:,2)<(rect(2)+rect(4)));
    idx_storm_include = 1:length(brp_all);
    storm_outlier_idx = cat(1,storm_outlier_idx,outlier_idx_now);
%     idx_storm_include = 1:length(brp_all); 
    idx_storm_include(intersect(idx_storm_include,storm_outlier_idx))=[];
    clf
    plot(brp_all(idx_storm_include,1),brp_all(idx_storm_include,2),'b+','MarkerSize',5)
    axis equal
end


%% compute hull
clf
dist_thresh = 60;
[hull_idx] = compute_hull(brp_all(idx_storm_include,:),dist_thresh);
clf
hold on

plot_hull(brp_all(idx_storm_include,:),hull_idx)

% plot(brp_all(idx_storm_include,1),brp_all(idx_storm_include,2),'m+')

% find distance between points and edges
hull_edge_1 = brp_all(idx_storm_include(hull_idx),:);
idx2 = [hull_idx(2:end) hull_idx(1)];
hull_edge_2 = brp_all(idx_storm_include(idx2),:);

brp_use = brp_all(idx_storm_include,:);

for i=1:length(brp_use)
    point = brp_use(i,:);
    for j = 1:length(hull_idx)
        line1 = hull_edge_1(j,:);
        line2 = hull_edge_2(j,:);
        dist_edge_mat(j) = point_to_line(point,line1,line2);
    end
    brp_dist(i) = min(dist_edge_mat);
end

 %plot interior points
clf
plot_hull(brp_use,hull_idx)
interior_idx = find(brp_dist>0);
hold on
 clf
[xxx,yyy] = points_to_spline(hull_idx,brp_use,0.95);
plot(brp_use(:,1),brp_use(:,2),'k.')
hold on
% plot(brp_use(interior_idx,1),brp_use(interior_idx,2),'m+')

%% save output
save('hull1','xxx','yyy')


%% quantify


for i=1:length(brp_use)
    idx_original = idx_storm_include(i); %original brp idx
    idx_match = find(out.match_idx(:,1)==idx_original); %brp match idx
    if ~~idx_match %check if match exists
        dist_match_brp(i) = brp_dist(i);
        heatmap_match_idx = out.match_idx(idx_match,2); %fcn match idx
        heatmap_count(i) = numel(find(out.fcn_cluster_idx==heatmap_match_idx));
    else
        dist_match_brp(i) = brp_dist(i);
        heatmap_count(i) = 0;
    end
end

clf
plot(dist_match_brp,heatmap_count,'k.')
%%
nomatch_idx = find(heatmap_count==0);
bins = 0:2:20;
[counts1 edge1 mid1 loc1] = histcn(dist_match_brp(nomatch_idx)',bins);
[counts2 edge2 mid2 loc2] = histcn(dist_match_brp',bins);
count_norm = counts1./counts2;
count_norm(isnan(count_norm))=0;
mid1=cell2mat(mid1);
% 
% std = std(total_norm,0,2)./sqrt(7);
plot(mid1(1:9),count_norm(1:9)
