%do affine transformation, clustering, and point matching

addpath ../common
addpath ../clustering
addpath ../image-segmentation
%load in data first
% x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
% y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];

%subsample - load heatmap then storm data
[r_hm,filehead]=OpenMolList;
[r,filehead]=OpenMolList;

brp_idx = find(r.cat==2);
cac_idx = find(r.cat==1);

%% get rectangle for cropping
clf
plot(r_hm.xc,r_hm.yc,'k.')
% ax=gca
% rect=getrect(ax)
% hold on
% heatmap_idx_crop = find(r.cat==1&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));
% plot(r.xc(heatmap_idx_crop),r.yc(heatmap_idx_crop),'m.')
% x1_heatmap = r.xc(heatmap_idx_crop);
% y1_heatmap = r.yc(heatmap_idx_crop);


% brp_idx = find(r.cat==3&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));
% 
% cac_idx = find(r.cat==2&r.xc>rect(1)&r.yc>rect(2)&...
%     r.xc<(rect(1)+rect(3))&r.yc<(rect(2)+rect(4)));

%% function [centers, counts] = dbscan_fcn(x,y,MinPts,eps);
%note - cluster #1 is non-clustered data
[centers_heatmap,counts_heatmap,cluster_idx_heatmap] = dbscan_fcn(r_hm.xc,r_hm.yc,15,1.5);
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
clf
plot(r_hm.xc,r_hm.yc,'k.')
axis equal
hold on
plot(centers_heatmap(1:end,1),centers_heatmap(1:end,2),'m+','MarkerSize',5)
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
    plot(r_hm.xc,r_hm.yc,'k.')
    axis equal
    hold on
    plot(centers_heatmap(2:end,1),centers_heatmap(2:end,2),'m+','MarkerSize',5)
    plot(x_brp,y_brp,'g.')
end



%% STORM peaks w/ HistTxt2D

crop_size = 3;
KeepGoing=1;
X=[x_brp y_brp];
% while KeepGoing==1
%     prompt = 'Enter peak threshold (press 0 to exit): ';
%     thresh = input(prompt);
%     if thresh==0
%         cpreg=0;
%         break
%     end
%%
    f=fspecial('disk',2);
    g = fspecial('gauss',3);
%     [coords,count_copy,count_coords] = HistTxt2D_fcn(X,thresh,crop_size);
    
    
    
    [count1 edges1 mid1 loc1] = histcn(X, 0:ceil(max(x_brp)), 0:ceil(max(y_brp)));
    clf
    count_copy_filter = imfilter(count1,f)./max(count1(:));
    imshow(imadjust(count_copy_filter))
    hold on
%     plot(count_coords(:,2),count_coords(:,1),'m+','MarkerSize',5)
    [row,col] = find(imregionalmax(count_copy_filter));
    plot(col,row,'m+','MarkerSize',5)
    
   %%

end
%% crop outliers
% keep track of included and deleted indices so matches can be referenced back to
% original list
KeepGoing=1;
idx_hm_include = 1:length(centers_heatmap);
idx_storm_include = 1:length(coords);
hm_outlier_idx = [];
storm_outlier_idx = [];


while KeepGoing==1
    prompt = ('Keep going? 1 to continue, 0 to exit')
    KeepGoing=input(prompt);
    if KeepGoing==0
        break
    end    
    clf
    plot(centers_heatmap(idx_hm_include,1),centers_heatmap(idx_hm_include,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(coords(idx_storm_include,1),coords(idx_storm_include,2),'b+','MarkerSize',5)
    ax=gca;
    rect=getrect(ax)
    hm_outlier_idx_now = find(centers_heatmap(:,1)>rect(1)&centers_heatmap(:,2)>rect(2)&...
    centers_heatmap(:,1)<(rect(1)+rect(3))&centers_heatmap(:,2)<(rect(2)+rect(4)));
    storm_outlier_idx_now = find(coords(:,1)>rect(1)&coords(:,2)>rect(2)&...
    coords(:,1)<(rect(1)+rect(3))&coords(:,2)<(rect(2)+rect(4)));
    hm_outlier_idx = cat(1,hm_outlier_idx,hm_outlier_idx_now);
    storm_outlier_idx = cat(1,storm_outlier_idx,storm_outlier_idx_now);
    % delete from list
    idx_hm_include = 1:length(centers_heatmap);
    idx_storm_include = 1:length(coords); 
    idx_hm_include(intersect(idx_hm_include,hm_outlier_idx))=[];
    idx_storm_include(intersect(idx_storm_include,storm_outlier_idx))=[];
    clf
    plot(centers_heatmap(idx_hm_include,1),centers_heatmap(idx_hm_include,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(coords(idx_storm_include,1),coords(idx_storm_include,2),'b+','MarkerSize',5)
end


%% do point matching on defined ROIs
storm_coords_include = coords(idx_storm_include,:);
heatmap_coords_include = centers_heatmap(idx_hm_include,:);

hm_delete = [];
storm_delete = [];
storm_match_final = [];
hm_match_final = [];

idx_hm_include_match = 1:length(heatmap_coords_include);
idx_storm_include_match = 1:length(storm_coords_include);

PointMatch=1;
while ~~PointMatch
    
    clf
    plot(heatmap_coords_include(idx_hm_include_match,1),...
        heatmap_coords_include(idx_hm_include_match,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(storm_coords_include(idx_storm_include_match,1),...
        storm_coords_include(idx_storm_include_match,2),'b+','MarkerSize',5)
    ax=gca;
    rect=getrect(ax)
    %find points within rectangle
    hm_match = find(heatmap_coords_include(idx_hm_include_match,1)>rect(1)&...
                heatmap_coords_include(idx_hm_include_match,2)>rect(2)&...
                heatmap_coords_include(idx_hm_include_match,1)<(rect(1)+rect(3))&...
                heatmap_coords_include(idx_hm_include_match,2)<(rect(2)+rect(4)));

    %also find points that haven't been matched yet
    
    hm_match_2 = find(~ismember(heatmap_coords_include,hm_delete));

    storm_match = find(storm_coords_include(idx_storm_include_match,1)>rect(1)&...
                storm_coords_include(idx_storm_include_match,2)>rect(2)&...
                storm_coords_include(idx_storm_include_match,1)<(rect(1)+rect(3))&...
                storm_coords_include(idx_storm_include_match,2)<(rect(2)+rect(4)));
            
%             clf
%     plot(heatmap_coords_include(hm_match,1),heatmap_coords_include(hm_match,2),'m+','MarkerSize',5)
%     hold on
%         plot(storm_coords_include(hm_match,1),storm_coords_include(hm_match,2),'m+','MarkerSize',5)
    
    %remove these indices for next step
    hm_delete = cat(2,hm_delete,idx_hm_include_match(hm_match));
    storm_delete = cat(2,storm_delete,idx_storm_include_match(storm_match));


    
    [match_idx_now] = pointMatch(storm_coords_include(idx_storm_include_match(storm_match),:),...
                    heatmap_coords_include(idx_hm_include_match(hm_match),:));
    %index back to include_ list
    storm_match_now = idx_storm_include_match(storm_match(match_idx_now(:,1)));
    hm_match_now = idx_hm_include_match(hm_match(match_idx_now(:,2)));
    
    storm_match_final = cat(2,storm_match_final,storm_match_now);
    hm_match_final = cat(2,hm_match_final,hm_match_now);
    
    idx_hm_include_match = 1:length(heatmap_coords_include);
    idx_storm_include_match = 1:length(storm_coords_include); 
    idx_hm_include_match(intersect(idx_hm_include_match,hm_delete))=[];
    idx_storm_include_match(intersect(idx_storm_include_match,storm_delete))=[];
    
    
%     keyboard
    
end
%% do point matching
% need to write temporary point files in path
% set parameter -t for local geometry preservation
% 
% for i=0.01:0.01:0.1


clf
plot(heatmap_coords_include(hm_match_final,1),heatmap_coords_include(hm_match_final,2),'m+','MarkerSize',5)
hold on
axis equal
plot(storm_coords_include(storm_match_final,1),storm_coords_include(storm_match_final,2),'b+','MarkerSize',5)


% clf
% plot(centers_heatmap(:,1),centers_heatmap(:,2),'m+','MarkerSize',5)
% hold on
% axis equal
% plot(coords(:,1),coords(:,2),'b+','MarkerSize',5)

match_idx = cat(1,storm_match_final,hm_match_final)';

idx_storm_match_original = idx_storm_include(match_idx(:,1));
idx_hm_match_original = idx_hm_include(match_idx(:,2));

outlier_thresh=1500;
for j=1:length(hm_match_final)
    %referenced back to originial uncropped list
    heatmap_idx_now = idx_hm_include(match_idx(j,2));
    storm_idx_now = idx_storm_include(match_idx(j,1));
    
    heatmap_idx_now = idx_hm_match_original(j);
    storm_idx_now = idx_storm_match_original(j);
    
    dist(j) = l2_dist_mat(centers_heatmap(heatmap_idx_now,:)',coords(storm_idx_now,:)');
    x=[centers_heatmap(heatmap_idx_now,1) coords(storm_idx_now,1)];
    y=[centers_heatmap(heatmap_idx_now,2) coords(storm_idx_now,2)];
    if dist(j)>outlier_thresh
        line(x,y,'Color','red')
    else
        line(x,y)
    end
end
plot_out=sprintf('C:\\STORM\\image_registration\\plots\\plot_%2.2f.png',i)
%     print(gcf,plot_out,'-dpng');
% end
%%


match_idx_filter = cat(1,idx_storm_match_original,idx_hm_match_original)';

idx_delete = find(dist>outlier_thresh);

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
%%
out.cac_count_match = cac_count(match_idx_filter(:,1));
out.brp_count_match = brp_count(match_idx_filter(:,1));
out.heatmap_count_match = counts_heatmap(match_idx_filter(:,2));
out.heatmap_count_match = out.heatmap_count_match';

cac_match_idx = cac_count.*0;
cac_match_idx(match_idx_filter(:,1))=1;

cac_nomatch_idx = find(cac_match_idx==0);
out.cac_nomatch = cac_count(cac_nomatch_idx);

out.brp_nomatch = brp_count(cac_nomatch_idx);
%%

coords_test = coords(match_idx_filter(:,1),:);
hm_test = centers_heatmap(match_idx_filter(:,2),:);
clf
plot(coords_test(:,1),coords_test(:,2),'b.')
hold on
plot(hm_test(:,1),hm_test(:,2),'m.')