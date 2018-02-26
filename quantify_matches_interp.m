%plot matches

%1. compute interpolation from cleaned matches
%2. compute hull from out struct
%load:
%---saved output from interpolate_matches (points_interp),
%---out struct from initial matching (for original warped heatmap coordinates)
%---hull (compute_mask)
addpath ../common
%%
storm_centers = out.storm_centers;

[r,filehead]=OpenMolList;

brp_idx = find(r.cat==2);
cac_idx = find(r.cat==1&r.valid>2000);       %added filtering
brp_cac_idx = cat(1,brp_idx,cac_idx);

brp_coords = [r.xc(brp_idx),r.yc(brp_idx)];
cac_coords = [r.xc(cac_idx),r.yc(cac_idx)];
fcn_coords = out.fcn_loc_coords;

%% quantify each cluster

brp_dist_mat = l2_dist(storm_centers,brp_coords);
cac_dist_mat = l2_dist(storm_centers,cac_coords);
fcn_dist_mat = l2_dist(fcn_center_predict,fcn_coords);

dist_thresh = 1.5;
dist_thresh_sq = dist_thresh.^2;

brp_num=[];
cac_num=[];
fcn_num=[];
for i=1:length(storm_centers)
    brp_num(i)=numel(find(brp_dist_mat(i,:)<dist_thresh_sq));
    cac_num(i)=numel(find(cac_dist_mat(i,:)<dist_thresh_sq));
end

for i=1:length(storm_centers)
    fcn_num(i)=numel(find(fcn_dist_mat(i,:)<dist_thresh_sq));
 
end

hullx = xxx;
hully = yyy;
edge_dist_mat = min(l2_dist(storm_centers,[xxx yyy]),[],2);


%% crop outliers
% idx_hm_include = 1:length(fcn_centers);
idx_storm_include = 1:length(storm_centers);

% fcn_centers_include_idx=1:length(match_idx);

hm_outlier_idx = [];
storm_outlier_idx = [];

select_rect=1;
while select_rect==1
    clf
    plot(fcn_center_predict(idx_storm_include,1),fcn_center_predict(idx_storm_include,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(xxx,yyy,'c.','LineWidth',2)
    plot(storm_centers(idx_storm_include,1),storm_centers(idx_storm_include,2),'b+','MarkerSize',5)
%     for j=1:length(match_idx_crop)
%     %referenced back to originial uncropped list
%         heatmap_idx_now = match_idx_crop(j,2);
%         storm_idx_now = match_idx_crop(j,1);
% 
%         x=[fcn_centers(heatmap_idx_now,1) storm_centers(storm_idx_now,1)];
%         y=[fcn_centers(heatmap_idx_now,2) storm_centers(storm_idx_now,2)];
% 
%         line(x,y)
% 
%     end
    ax=gca;
    
    try
        rect=getrect(ax);
    catch
        select_rect=0;
        break
    end

%     hm_outlier_idx_now = find(fcn_centers(:,1)>rect(1)&fcn_centers(:,2)>rect(2)&...
%     fcn_centers(:,1)<(rect(1)+rect(3))&fcn_centers(:,2)<(rect(2)+rect(4)));
    storm_centers_outliers_temp = find(storm_centers(:,1)>rect(1)&...
                                storm_centers(:,2)>rect(2)&...
                                storm_centers(:,1)<(rect(1)+rect(3))&...
                                storm_centers(:,2)<(rect(2)+rect(4)));
    idx_storm_include = 1:length(storm_centers);                         
    storm_outlier_idx = cat(1,storm_outlier_idx,storm_centers_outliers_temp);
    storm_outlier_ids = unique(storm_outlier_idx);
    idx_storm_include(storm_outlier_idx)=[];
    
end

%% quantify distance dependence
fcn_num_use = fcn_num(idx_storm_include);
brp_num_use = brp_num(idx_storm_include);
cac_num_use = cac_num(idx_storm_include);
edge_dist_mat_sqrt = sqrt(edge_dist_mat(idx_storm_include));
clf
abc = fcn_num_use./brp_num_use;
plot(edge_dist_mat_sqrt,fcn_num_use,'m.')

out_mat = [idx_storm_include' edge_dist_mat_sqrt...
            fcn_num_use' brp_num_use' abc'];

        
num_bins = 4;
pr=fcn_num_use/200;

% [pr_sort pr_sort_ind]=sort(pr);
[d_sort d_sort_ind]=sort(edge_dist_mat_sqrt);
pr_sort = pr(d_sort_ind);
brp_sort = brp_num_use(d_sort_ind);
cac_sort = cac_num_use(d_sort_ind);
bin_num=round(length(pr_sort)/num_bins);

n=sqrt(length(idx_storm_include));
for i = 1:num_bins
    lim1 = bin_num*(i-1)+1;
    lim2 = min(bin_num*i,numel(brp_sort));
    pr_binned(i) = mean(pr_sort(lim1:lim2));
    pr_std(i) = std(pr_sort(lim1:lim2))/n;
    dist_binned(i) = mean(d_sort(lim1:lim2));
    dist_std(i) = std(d_sort(lim1:lim2))/n;
    brp_binned(i) = mean(brp_sort(lim1:lim2));
    brp_std(i) = std(brp_sort(lim1:lim2))/n;
    cac_binned(i) = mean(cac_sort(lim1:lim2));
    cac_std(i) = std(cac_sort(lim1:lim2))/n;

end     

clf

out_mat_2 = [pr_binned' pr_std' brp_binned' brp_std' cac_binned' cac_std' dist_binned' dist_std'];

% plot(pr_binned,cac_binned,'m+')
MarkerSize = 5;
   clf
        hold on
        for q=1:num_bins
            plot(dist_binned,pr_binned,'m.')
        end
        errorbarxy(dist_binned,pr_binned,dist_std,pr_std,[1 0 0]);
                                
                            


%         set(leg,'location','NorthWest')
        XLimits=xlim;
        YLimits=ylim;
        xlim([0,max(dist_binned)*1.2])
        ylim([0,YLimits(2)])
        xlabel('Distance to edge')
        ylabel('Pr')
        FigureStandardizer(get(gca,'xlabel'), get(gca,'ylabel'), gca);

        set(gcf, 'color', 'white');
%         set(leg,'location','EastOutside')
        set(gcf,'units','normalized','position',[0.1,0.1,0.6,0.8])
        set(gca,'position',[0.15,0.25,0.5,0.5])
% axis equal
% figure
% plot(cac_num,fcn_num,'m.')
% axis equal
% 
%     plot(fcn_centers(:,1),fcn_centers(:,2),'m.','MarkerSize',3)
%     hold on
%     axis equal
%     plot(storm_centers(:,1),storm_centers(:,2),'b.','MarkerSize',3)
%     plot(storm_centers(storm_centers_include_idx,1),storm_centers(storm_centers_include_idx,2),'k+')
%     hold on
%     plot(fcn_centers(fcn_centers_include_idx,1),fcn_centers(fcn_centers_include_idx,2),'m+')

%% plot and save
    figure
    plot(fcn_centers(match_idx(:,2),1),fcn_centers(match_idx(:,2),2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(storm_centers(match_idx(:,1),1),storm_centers(match_idx(:,1),2),'b+','MarkerSize',5)
    for j=1:length(match_idx)
    %referenced back to originial uncropped list
        heatmap_idx_now = match_idx(j,2);
        storm_idx_now = match_idx(j,1);

        x=[fcn_centers(heatmap_idx_now,1) storm_centers(storm_idx_now,1)];
        y=[fcn_centers(heatmap_idx_now,2) storm_centers(storm_idx_now,2)];

        line(x,y,'Color','r')

    end

    plot(fcn_centers(match_idx_crop(:,2),1),fcn_centers(match_idx_crop(:,2),2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(storm_centers(match_idx_crop(:,1),1),storm_centers(match_idx_crop(:,1),2),'b+','MarkerSize',5)
    for j=1:length(match_idx_crop)
    %referenced back to originial uncropped list
        heatmap_idx_now = match_idx_crop(j,2);
        storm_idx_now = match_idx_crop(j,1);

        x=[fcn_centers(heatmap_idx_now,1) storm_centers(storm_idx_now,1)];
        y=[fcn_centers(heatmap_idx_now,2) storm_centers(storm_idx_now,2)];

        line(x,y,'Color','b')

    end
    
filename_plot = sprintf('%s_match.png',filehead);
print(gcf,filename_plot,'-dpng','-r600');

filename_match = sprintf('%s_match_idx',filehead);
save(filename_match,'match_idx_crop')
    
    