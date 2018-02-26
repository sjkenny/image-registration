% compute density of STORM clusters in neuron


addpath ../common
%% if data needs to be clustered
[r,filename]=OpenMolList;   

brp_cat = 2;
brp_idx = find(r.cat==brp_cat);


PxSize=0.2;
scale=1/PxSize;
centers_heatmap_scale=centers_heatmap.*scale;

KeepGoing=1;
while KeepGoing==1
    prompt = 'Enter density threshold (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
    brp_idx_filter = find(r.valid(brp_idx)>thresh);
    
    x_brp_filter=r.xc(brp_idx(brp_idx_filter));
    y_brp_filter=r.yc(brp_idx(brp_idx_filter));
    X = [y_brp_filter x_brp_filter];
    [count1 edges1 mid1 loc1] = histcn(X, 0:PxSize:ceil(max(y_brp_filter)), 0:PxSize:ceil(max(x_brp_filter)));
    
%     [im_ff,center_list] = flood_fill_indices(count1);
    f=fspecial('disk',0.5/PxSize);
    g = fspecial('gauss',3/PxSize);
    clf
    count_copy_filter = imfilter(count1,f);
    count_copy_filter2 = uint16(10000*(count_copy_filter./max(count_copy_filter(:))));
    
    [im_ff,storm_centers] = flood_fill_indices(count_copy_filter);
    
    imshow(imadjust(count_copy_filter))
    axis equal
    hold on
%     [row,col] = find(imregionalmax(count_copy_filter));
    plot(storm_centers(:,1),storm_centers(:,2),'b+','MarkerSize',5)

   
end
filename_out = sprintf('%s-brp-centers',filename);
save(filename_out,'storm_centers')
%% 



