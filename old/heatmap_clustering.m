%do clustering stuff on heatmap data
%define ROI of 3-color data with heatmap and save .txt file

addpath ../common
addpath ../clustering
addpath ../image-segmentation
%load in data first
% x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
% y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];

%subsample
[r,filehead]=OpenMolListTxt;


heatmap_idx = find(r.cat==1);
brp_idx = find(r.cat==3);
cac_idx = find(r.cat==2);

%% get rectangle for cropping
clf
plot(r.xc(heatmap_idx),r.yc(heatmap_idx),'k.')
% ax=gca
% rect=getrect(ax)
% hold on
% heatmap_idx_crop = find(r.cat==1&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));
% plot(r.xc(heatmap_idx_crop),r.yc(heatmap_idx_crop),'m.')
% x1_heatmap = r.xc(heatmap_idx_crop);
% y1_heatmap = r.yc(heatmap_idx_crop);
x1_heatmap = r.xc(heatmap_idx);
y1_heatmap = r.yc(heatmap_idx);

% brp_idx = find(r.cat==3&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));
% 
% cac_idx = find(r.cat==2&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));

%% function [centers, counts] = dbscan_fcn(x,y,MinPts,eps);
%note - cluster #1 is non-clustered data
[centers_heatmap,counts_heatmap,cluster_idx_heatmap] = dbscan_fcn(x1_heatmap,y1_heatmap,8,1.2);
%delete centers_heatmap(1)
centers_heatmap(1,:)=[];
counts_heatmap(1)=[];
cluster_idx_heatmap=cluster_idx_heatmap-1;
% for i = 2:length(centers)
%     
%     clf
%     plot(x1,y1,'k.')
%     axis equal
%     hold on
%     idx_use = find(idx==i);
%     plot(x1(idx_use),y1(idx_use),'m.','MarkerSize',15)
%     plot(centers(2:end,1),centers(2:end,2),'m+','MarkerSize',5)
%     plot(centers(i,1),centers(i,2),'b+','MarkerSize',15)
%     keyboard
% end
%%
clf
plot(x1_heatmap,y1_heatmap,'k.')
axis equal
hold on
plot(centers_heatmap(2:end,1),centers_heatmap(2:end,2),'m+','MarkerSize',5)
sprintf('%d',length(centers_heatmap)-1)
%%
% clf
% % imshow(imadjust(Probability_Heatmap_Filtered))
% hold on
% plot(x1,y1,'k.')
% plot(centers_heatmap(2:end,1),centers_heatmap(2:end,2),'m+','MarkerSize',10)


%% storm filtering

% r=OpenMolListTxt;
% 
% 
KeepGoing=1;
while KeepGoing==1
    prompt = 'Enter density threshold (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
    brp_idx_filter = find(r.valid(brp_idx)>thresh);
    
    x_brp=r.xc(brp_idx(brp_idx_filter));
    y_brp=r.yc(brp_idx(brp_idx_filter));
    clf
    plot(x1_heatmap,y1_heatmap,'k.')
    axis equal
    hold on
    plot(centers_heatmap(2:end,1),centers_heatmap(2:end,2),'m+','MarkerSize',5)
    plot(x_brp,y_brp,'g.')
end



%% STORM peaks w/ HistTxt2D

crop_size = 3;
KeepGoing=1;
X=[x_brp y_brp];
while KeepGoing==1
    prompt = 'Enter peak threshold (press 0 to exit): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
    [coords,count_copy,count_coords] = HistTxt2D(X,thresh,crop_size);
    clf
    imshow(count_copy)
    hold on
    plot(count_coords(:,2),count_coords(:,1),'m+','MarkerSize',5)

end


%% do point matching
% need to write temporary point files in path
% set parameter -t for local geometry preservation

[outfile_storm outfile_heatmap] = write_Txt(filehead,coords,centers_heatmap);
exe_name = 'PointMatchDemo2.exe';
system(sprintf('%s -t 0.2 -f -i 5 storm_coords.point heatmap_coords.point',exe_name))

% match centers
A=importdata('storm_coords_heatmap_coords.match');
match_idx = A.data;
match_idx=match_idx+1;
clf
plot(centers_heatmap(:,1),centers_heatmap(:,2),'m+','MarkerSize',5)
hold on
plot(coords(:,1),coords(:,2),'b+','MarkerSize',5)
%%
outlier_thresh=400;
for i=1:length(match_idx)
    heatmap_idx_now = match_idx(i,2);
    storm_idx_now = match_idx(i,1);
    dist(i) = l2_dist_mat(centers_heatmap(heatmap_idx_now,:)',coords(storm_idx_now,:)');
    x=[centers_heatmap(heatmap_idx_now,1) coords(storm_idx_now,1)];
    y=[centers_heatmap(heatmap_idx_now,2) coords(storm_idx_now,2)];
    if dist(i)>outlier_thresh
        line(x,y,'Color','red')
    else
        line(x,y)
    end
end
idx_delete = find(dist>outlier_thresh);
match_idx_filter = match_idx;
match_idx_filter(idx_delete,:)=[];

% get statistics from matched locs
% ratio of Cac within each cluster to # of heatmap locs
% need: # of heatmap locs per cluster (counts_heatmap)
%      # of cac locs within each cluster - normalize?
%       matched vs. non-matched Cac


%% count cac in each cluster
search_radius = 1;
cac_x = r.xc(cac_idx);
cac_y = r.yc(cac_idx);
brp_x = r.xc(brp_idx);
brp_y = r.yc(brp_idx);
dist_mat_cac = l2_dist_mat(coords',[cac_x cac_y]');
dist_mat_brp = l2_dist_mat(coords',[brp_x brp_y]');
cac_count=zeros(length(coords),1);
brp_count=zeros(length(coords),1);
for i=1:length(coords)
    cac_count(i) = numel(find(dist_mat_cac(:,i)<search_radius^2));
    brp_count(i) = numel(find(dist_mat_brp(:,i)<search_radius^2));
end

cac_count_match = cac_count(match_idx_filter(:,1));
brp_count_match = brp_count(match_idx_filter(:,1));
heatmap_count_match = counts_heatmap(match_idx_filter(:,2));
heatmap_count_match = heatmap_count_match';

cac_match_idx = cac_count.*0;
cac_match_idx(match_idx_filter(:,1))=1;

cac_nomatch_idx = find(cac_match_idx==0);
cac_nomatch = cac_count(cac_nomatch_idx);

brp_nomatch = brp_count(cac_nomatch_idx);

