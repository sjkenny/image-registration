% compute density of STORM clusters in neuron


addpath ../common
%% if data needs to be clustered
[r,filename]=OpenMolList;   
%% density based thresholding and flood fill clustering
brp_cat = 2;
brp_idx = find(r.cat==brp_cat);

scale = 5; %upscale binned image for flood fill
PxSize = 1/scale;


KeepGoing=1;
while KeepGoing==1
    prompt = 'Enter density threshold (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
    % get XY coords of STORM image for binning
    brp_idx_filter = find(r.valid(brp_idx)>thresh); %minimum local density threshold
    x_brp_filter=r.xc(brp_idx(brp_idx_filter));
    y_brp_filter=r.yc(brp_idx(brp_idx_filter));
    X = [y_brp_filter x_brp_filter];
    %2D binning 
    [count1 edges1 mid1 loc1] = histcn(X, 0:PxSize:ceil(max(y_brp_filter)), 0:PxSize:ceil(max(x_brp_filter)));
    
    %2 step filtering: disk and then gaussian for cluster amplification
    f=fspecial('disk',0.5/PxSize);
    g = fspecial('gauss',2/PxSize);
    clf
    i1 = imfilter(count1,f);
    i2 = imfilter(i1,g);
    i3 = i2./max(i2(:));
    imgbw = im2bw(i3,1.8*mean(i3(i3>0)));
   
    [im_ff,center_list] = flood_fill_indices(imgbw);
    imshow(imgbw)
    [im_ff,storm_centers_scale,area_idx] = flood_fill_indices(imgbw);
    imgbw_perim = bwperim(imgbw);
    
%     imshow(imgbw_perim) %for STORM overlay
    axis equal
    hold on
%     [row,col] = find(imregionalmax(count_copy_filter)); gives multiple
%     maxima per cluster
    plot(storm_centers_scale(:,1),storm_centers_scale(:,2),'b+','MarkerSize',5)

   
end

storm_centers=storm_centers_scale.*PxSize;  %rescale to original coordinates
filename_out = sprintf('%s-brp-centers',filename);
save(filename_out,'storm_centers')


for i=1:length(area_idx)
    area_list(i) = length(area_idx{i});
end
%%
area_list_norm = area_list.*PxSize^2;
area_bins = 0:15;
area_bins = exp(0.3*area_bins)-1;
[count_area edges_area mid_area loc_area] = histcn(area_list_norm',area_bins');
plot(cell2mat(mid_area),count_area,'k')
%% crop outliers
idx_include = crop_outliers(storm_centers);

%% compute hull

storm_centers_include = storm_centers(idx_include,:);
clf
dist_thresh = 70;
[hull_idx] = compute_hull(storm_centers_include,dist_thresh);
clf
hold on

plot_hull(storm_centers_include,hull_idx)

hull_idx_wrap = cat(2,hull_idx, hull_idx(1));   %wrap onto first point

hull_area = polyarea(storm_centers_include(hull_idx_wrap,1),...
                     storm_centers_include(hull_idx_wrap,2));
                 
global_density = length(storm_centers_include)/hull_area;

