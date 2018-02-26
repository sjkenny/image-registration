%plot matches
%load in output
addpath ../common
%%
fcn_centers = out.fcn_centers;
storm_centers = out.storm_centers;
match_idx = out.match_idx;

% plot(fcn_centers(:,1),fcn_centers(:,2),'m.','MarkerSize',3)
% hold on
% axis equal
% plot(storm_centers(:,1),storm_centers(:,2),'b.','MarkerSize',3)

[r,filehead]=OpenMolList;

brp_idx = find(r.cat==2);
cac_idx = find(r.cat==1&r.valid>2000);       %added filtering
brp_cac_idx = cat(1,brp_idx,cac_idx);

brp_coords = [r.xc(brp_idx),r.yc(brp_idx)];
cac_coords = [r.xc(cac_idx),r.yc(cac_idx)];
fcn_coords = out.fcn_loc_coords;

%% crop outliers
idx_hm_include = 1:length(fcn_centers);
idx_storm_include = 1:length(storm_centers);

fcn_centers_include_idx=1:length(match_idx);

hm_outlier_idx = [];
storm_outlier_idx = [];
match_idx_crop = match_idx;
select_rect=1;
while select_rect==1
    clf
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

        line(x,y)

    end
    ax=gca;
    
    try
        rect=getrect(ax);
    catch
        select_rect=0;
        break
    end
    match_idx_crop=match_idx;
%     hm_outlier_idx_now = find(fcn_centers(:,1)>rect(1)&fcn_centers(:,2)>rect(2)&...
%     fcn_centers(:,1)<(rect(1)+rect(3))&fcn_centers(:,2)<(rect(2)+rect(4)));
    storm_centers_outliers_temp = find(storm_centers(:,1)>rect(1)&...
                                storm_centers(:,2)>rect(2)&...
                                storm_centers(:,1)<(rect(1)+rect(3))&...
                                storm_centers(:,2)<(rect(2)+rect(4)));
                            
    a=find(ismember(match_idx_crop(:,1),storm_centers_outliers_temp));    
    storm_outlier_idx = cat(1,storm_outlier_idx,a);
    match_idx_crop(storm_outlier_idx,:)=[];
    
end
%% save cleaned matches
filename_match = sprintf('%s_match_idx',filehead);
save(filename_match,'match_idx_crop')

%% select regions to quantify
% storm_centers_brp_include = centers_brp(idx_storm_include,:);
% heatmap_centers_brp_include = centers_heatmap(idx_hm_include,:);

select_rect=1;
while ~~select_rect
    
    clf
    plot(fcn_centers(match_idx_crop(:,2),1),fcn_centers(match_idx_crop(:,2),2),'m.','MarkerSize',3)
    hold on
    axis equal
    plot(storm_centers(match_idx_crop(:,1),1),storm_centers(match_idx_crop(:,1),2),'b.','MarkerSize',3)
    for j=1:length(match_idx_crop)
    %referenced back to originial uncropped list
        heatmap_idx_now = match_idx_crop(j,2);
        storm_idx_now = match_idx_crop(j,1);

        x=[fcn_centers(heatmap_idx_now,1) storm_centers(storm_idx_now,1)];
        y=[fcn_centers(heatmap_idx_now,2) storm_centers(storm_idx_now,2)];

        line(x,y)
        ax=gca;
    end
    try
        rect=getrect(ax)
    catch
        select_rect=0;
        break
    
    end
    
    %find points within rectangle
    storm_centers_include_temp = find(storm_centers(match_idx_crop(:,1),1)>rect(1)&...
                                storm_centers(match_idx_crop(:,1),2)>rect(2)&...
                                storm_centers(match_idx_crop(:,1),1)<(rect(1)+rect(3))&...
                                storm_centers(match_idx_crop(:,1),2)<(rect(2)+rect(4)));
                            
%     a=find(ismember(match_idx_crop(:,1),storm_centers_include_temp));
    fcn_centers_include_idx = match_idx_crop(storm_centers_include_temp,2);
    storm_centers_include_idx = match_idx_crop(storm_centers_include_temp,1);
    


end



clf

storm_centers_include = storm_centers(storm_centers_include_idx,:);
fcn_centers_include = fcn_centers(fcn_centers_include_idx,:);

% plot(storm_centers_include(:,1),storm_centers_include(:,2),'b+')
% hold on
% axis equal
% plot(fcn_centers_include(:,1),fcn_centers_include(:,2),'m+')

brp_dist_mat = l2_dist_mat(storm_centers_include',brp_coords');
cac_dist_mat = l2_dist_mat(storm_centers_include',cac_coords');
fcn_dist_mat = l2_dist_mat(fcn_centers_include',fcn_coords');

dist_thresh = 2;
brp_num=[];
cac_num=[];
fcn_num=[];
for i=1:length(storm_centers_include)
    brp_num(i)=numel(find(brp_dist_mat(:,i)<dist_thresh));
    cac_num(i)=numel(find(cac_dist_mat(:,i)<dist_thresh));
end

for i=1:length(fcn_centers_include)
    fcn_num(i)=numel(find(fcn_dist_mat(:,i)<dist_thresh));
 
end



%   convert to quartile and bin
num_bins = 4;
pr=fcn_num/200;
[pr_sort pr_sort_ind]=sort(pr);
brp_sort = brp_num(pr_sort_ind);
cac_sort = cac_num(pr_sort_ind);
bin_num=round(length(pr_sort)/num_bins);
n=sqrt(length(storm_centers_include));
for i = 1:num_bins
    lim1 = bin_num*(i-1)+1;
    lim2 = min(bin_num*i,numel(brp_sort));
    pr_binned(i) = mean(pr_sort(lim1:lim2));
    pr_std(i) = std(pr_sort(lim1:lim2))/n;
    brp_binned(i) = mean(brp_sort(lim1:lim2));
    brp_std(i) = std(brp_sort(lim1:lim2))/n;
    cac_binned(i) = mean(cac_sort(lim1:lim2));
    cac_std(i) = std(cac_sort(lim1:lim2))/n;
end

clf

out_mat = [pr_binned' pr_std' brp_binned' brp_std' cac_binned' cac_std'];
%
%% plot(pr_binned,cac_binned,'m+')
MarkerSize = 5;
   clf
        hold on
        for q=1:num_bins
            plot(pr_binned,brp_binned,'m.')
        end
        errorbarxy(pr_binned,brp_binned,pr_std,brp_std,[1 0 0]);
                                
                            


%         set(leg,'location','NorthWest')
        XLimits=xlim;
        YLimits=ylim;
        xlim([0,1])
        ylim([0,YLimits(2)])
        xlabel('Pr')
        ylabel('#Brp Loc')
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
    
    