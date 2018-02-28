%do clustering and point matching
%this is for new data with different classes of functional localizations
%data should be registered and saved in previous step img_reg_similarity
%currently 4 localization classes:
%reg_evoked_Ib
%reg_evoked_Is
%reg_spont_Ib
%reg_spont_Is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 1: load data
%Step 2: cluster functional data (dbscan)
%Step 3: cluster STORM data (flood fill or dbscan)
%Step 4: manually crop clusters outside neuron
%Step 5: automated point matching
%Step 6: clean matches and save output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Step 1: load data
% -STORM data as .bin file
% -functional data as affine registered .mat file
addpath ../common
addpath ../clustering
addpath ../image-segmentation


%parse input



%open storm list
[r,filehead]=OpenMolList;

% single_list=0;
% try
%     [r_hm,filehead_hm]=OpenMolList;
%     fprintf('Double molecule lists loaded\r')
% catch
%     fprintf('Single molecule list loaded\r')
%     single_list=1;
% end
% 
% if ~single_list
    brp_idx = find(r.cat==2);
    cac_idx = find(r.cat==1);
    brp_cac_idx = cat(1,brp_idx,cac_idx);
%     idx_1b = find(r_hm.cat==1);
%     hm_x_use = r_hm.xc(idx_1b);               %if functional data is                                          
%     hm_y_use = r_hm.yc(idx_1b);               %included in .bin file 
% else
%     brp_idx = find(r.cat==3);
%     cac_idx = find(r.cat==2);
%     brp_cac_idx = cat(1,brp_idx,cac_idx);
%     %cat1 = 1s
%     hm_idx=find(r.cat==1);
%     hm_x_use = r.xc(hm_idx);
%     hm_y_use = r.yc(hm_idx);

hm_x_use = bx1;
hm_y_use = by1;

% for classification
% reg_evoked_Ib = [bx1(coords_idx_crop==1) by1(coords_idx_crop==1)];
% reg_evoked_Is = [bx1(coords_idx_crop==2) by1(coords_idx_crop==2)];
% reg_spont_Ib = [bx1(coords_idx_crop==3) by1(coords_idx_crop==3)];
% reg_spont_Is = [bx1(coords_idx_crop==4) by1(coords_idx_crop==4)];


%% Step 2: cluster functional data (dbscan)


[centers_heatmap,counts_heatmap,cluster_idx_heatmap] = dbscan_fcn(hm_x_use,hm_y_use,7,1.5);
%delete centers_heatmap(1)
centers_heatmap(1,:)=[];
counts_heatmap(1)=[];
cluster_idx_heatmap=cluster_idx_heatmap-1;

clf
plot(hm_x_use,hm_y_use,'k.')
axis equal
hold on
plot(centers_heatmap(:,1),centers_heatmap(:,2),'m+','MarkerSize',5)
sprintf('%d',length(centers_heatmap)-1)



%% Step 3: cluster STORM data (flood fill)


% r=OpenMolListTxt;
PxSize=0.2;
scale=1/PxSize;
centers_heatmap_scale=centers_heatmap.*scale;
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
    
    x_brp_filter=r.xc(brp_idx(brp_idx_filter));
    y_brp_filter=r.yc(brp_idx(brp_idx_filter));
    X = [y_brp_filter x_brp_filter];
    [count1 edges1 mid1 loc1] = histcn(X, 0:PxSize:ceil(max(y_brp_filter)), 0:PxSize:ceil(max(x_brp_filter)));
    
%     [im_ff,center_list] = flood_fill_indices(count1);
    f=fspecial('disk',0.5/PxSize);
    g = fspecial('gauss',2/PxSize);
    clf
    i1 = imfilter(count1,f);
    i2 = imfilter(i1,g);
    i3 = i2./max(i2(:));
    imgbw = im2bw(i3,1.8*mean(i3(i3>0)));
   
    [im_ff,center_list] = flood_fill_indices(imgbw);
    
    imshow(imgbw)
    axis equal
    hold on
%     [row,col] = find(imregionalmax(count_copy_filter));
    plot(center_list(:,1),center_list(:,2),'b+','MarkerSize',5)
    plot(centers_heatmap_scale(2:end,1),centers_heatmap_scale(2:end,2),'m+','MarkerSize',5)
   
end
%% Step 3: cluster STORM data (dbscan)

% [centers_brp,counts_brp,cluster_idx_brp] = dbscan_fcn(x_brp_filter,y_brp_filter,50,0.4);
% %delete centers_heatmap(1)
% centers_brp(1,:)=[];
% counts_brp(1)=[];
% cluster_idx_brp=cluster_idx_brp-1;
% imshow(imadjust(count_copy_filter))
% axis equal
% hold on
% plot(centers_brp(:,1),centers_brp(:,2),'b+','MarkerSize',5)

%%

 
%% Step 4: manually crop clusters outside neuron
% keep track of included and deleted indices so matches can be referenced back to
% original list

centers_brp=center_list.*PxSize;

KeepGoing=1;
idx_hm_include = 1:length(centers_heatmap);
idx_storm_include = 1:length(centers_brp);
hm_outlier_idx = [];
storm_outlier_idx = [];


while KeepGoing==1
    clf
    plot(centers_heatmap(idx_hm_include,1),centers_heatmap(idx_hm_include,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(centers_brp(idx_storm_include,1),centers_brp(idx_storm_include,2),'b+','MarkerSize',5)
    ax=gca;
    try
        rect=getrect(ax);
    catch
        KeepGoing=0;
        break
    end
    
    hm_outlier_idx_now = find(centers_heatmap(:,1)>rect(1)&centers_heatmap(:,2)>rect(2)&...
    centers_heatmap(:,1)<(rect(1)+rect(3))&centers_heatmap(:,2)<(rect(2)+rect(4)));
    storm_outlier_idx_now = find(centers_brp(:,1)>rect(1)&centers_brp(:,2)>rect(2)&...
    centers_brp(:,1)<(rect(1)+rect(3))&centers_brp(:,2)<(rect(2)+rect(4)));
    hm_outlier_idx = cat(1,hm_outlier_idx,hm_outlier_idx_now);
    storm_outlier_idx = cat(1,storm_outlier_idx,storm_outlier_idx_now);
    % delete from list
    idx_hm_include = 1:length(centers_heatmap);
    idx_storm_include = 1:length(centers_brp); 
    idx_hm_include(intersect(idx_hm_include,hm_outlier_idx))=[];
    idx_storm_include(intersect(idx_storm_include,storm_outlier_idx))=[];
    clf
    plot(centers_heatmap(idx_hm_include,1),centers_heatmap(idx_hm_include,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(centers_brp(idx_storm_include,1),centers_brp(idx_storm_include,2),'b+','MarkerSize',5)
end


%% Step 5: automated point matching on defined ROIs
storm_centers_brp_include = centers_brp(idx_storm_include,:);
heatmap_centers_brp_include = centers_heatmap(idx_hm_include,:);

hm_delete = [];
storm_delete = [];
storm_match_final = [];
hm_match_final = [];

idx_hm_include_match = 1:length(heatmap_centers_brp_include);
idx_storm_include_match = 1:length(storm_centers_brp_include);

PointMatch=1;
while ~~PointMatch
    
    clf
    plot(heatmap_centers_brp_include(idx_hm_include_match,1),...
        heatmap_centers_brp_include(idx_hm_include_match,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(storm_centers_brp_include(idx_storm_include_match,1),...
        storm_centers_brp_include(idx_storm_include_match,2),'b+','MarkerSize',5)
    ax=gca;
    
    try
        rect=getrect(ax)
    catch
        PointMatch=0;
        break
    
    end
    
    %find points within rectangle
    hm_match = find(heatmap_centers_brp_include(idx_hm_include_match,1)>rect(1)&...
                heatmap_centers_brp_include(idx_hm_include_match,2)>rect(2)&...
                heatmap_centers_brp_include(idx_hm_include_match,1)<(rect(1)+rect(3))&...
                heatmap_centers_brp_include(idx_hm_include_match,2)<(rect(2)+rect(4)));

    %also find points that haven't been matched yet
    
    hm_match_2 = find(~ismember(heatmap_centers_brp_include,hm_delete));

    storm_match = find(storm_centers_brp_include(idx_storm_include_match,1)>rect(1)&...
                storm_centers_brp_include(idx_storm_include_match,2)>rect(2)&...
                storm_centers_brp_include(idx_storm_include_match,1)<(rect(1)+rect(3))&...
                storm_centers_brp_include(idx_storm_include_match,2)<(rect(2)+rect(4)));
            
%             clf
%     plot(heatmap_centers_brp_include(hm_match,1),heatmap_centers_brp_include(hm_match,2),'m+','MarkerSize',5)
%     hold on
%         plot(storm_centers_brp_include(hm_match,1),storm_centers_brp_include(hm_match,2),'m+','MarkerSize',5)
    
    %remove these indices for next step
    hm_delete = cat(2,hm_delete,idx_hm_include_match(hm_match));
    storm_delete = cat(2,storm_delete,idx_storm_include_match(storm_match));


    
    [match_idx_now] = pointMatch(storm_centers_brp_include(idx_storm_include_match(storm_match),:),...
                    heatmap_centers_brp_include(idx_hm_include_match(hm_match),:));
    %index back to include_ list
    storm_match_now = idx_storm_include_match(storm_match(match_idx_now(:,1)));
    hm_match_now = idx_hm_include_match(hm_match(match_idx_now(:,2)));
    
    storm_match_final = cat(2,storm_match_final,storm_match_now);
    hm_match_final = cat(2,hm_match_final,hm_match_now);
    
    idx_hm_include_match = 1:length(heatmap_centers_brp_include);
    idx_storm_include_match = 1:length(storm_centers_brp_include); 
    idx_hm_include_match(intersect(idx_hm_include_match,hm_delete))=[];
    idx_storm_include_match(intersect(idx_storm_include_match,storm_delete))=[];
    
    
%     keyboard
    
end
%  plot matches

clf
% imshow(imadjust(count_copy_filter))

% hold on
% axis equal
% plot(hm_x_use,hm_y_use,'m.')
plot(centers_heatmap(:,1),centers_heatmap(:,2),'m.','MarkerSize',3)
hold on
axis equal
plot(centers_brp(:,1),centers_brp(:,2),'b.','MarkerSize',3)

plot(heatmap_centers_brp_include(hm_match_final,1),heatmap_centers_brp_include(hm_match_final,2),'m+','MarkerSize',3)
plot(storm_centers_brp_include(storm_match_final,1),storm_centers_brp_include(storm_match_final,2),'b+','MarkerSize',3)

match_idx = cat(1,storm_match_final,hm_match_final)';

idx_storm_match_original = idx_storm_include(match_idx(:,1));
idx_hm_match_original = idx_hm_include(match_idx(:,2));

% outfile_png = sprintf('%s-match_all',filehead_hm);
% print(outfile_png,gcf,'-dpng')

%distance threshold for outliers
outlier_thresh=1500;
%draw lines b/w matches
dist_mat=[];
for j=1:length(hm_match_final)
    %referenced back to originial uncropped list
    heatmap_idx_now = idx_hm_match_original(j);
    storm_idx_now = idx_storm_match_original(j);
    
    dist_mat(j) = l2_dist_mat(centers_heatmap(heatmap_idx_now,:)',centers_brp(storm_idx_now,:)');
    x=[centers_heatmap(heatmap_idx_now,1) centers_brp(storm_idx_now,1)];
    y=[centers_heatmap(heatmap_idx_now,2) centers_brp(storm_idx_now,2)];
    if dist_mat(j)>outlier_thresh
        line(x,y,'Color','red')
    else
        line(x,y)
    end
end
% filename_plot = sprintf('%s_match.png',filehead_hm);
% print(gcf,filename_plot,'-dpng','-r600');
% end
%% filter outliers
match_idx_use = cat(1,idx_storm_match_original,idx_hm_match_original)';
match_idx_filter = clean_matches(centers_heatmap(match_idx_use(:,2),:),centers_brp(match_idx_use(:,1),:));
match_idx_out = match_idx_use(match_idx_filter,:);

%% save output matches

[x1,y1]=tforminv(tform,centers_heatmap(:,1),centers_heatmap(:,2));  %warp back to original axes

out.fcn_loc_coords = [bx1 by1]; %registered heatmap coordinates

out.fcn_centers = centers_heatmap;
out.fcn_centers_original = [x1 y1];
out.fcn_class = coords_idx_crop; %Ib/Is/evoked/spont
out.fcn_cluster_idx = cluster_idx_heatmap;
out.storm_centers = centers_brp;
out.match_idx = match_idx_out;
out.tform = tform;
out.idx_storm_clean = idx_storm_include;    % match_idx references these
out.idx_fcn_clean = idx_hm_include;



filenameStats = sprintf('%s-stats1.mat',filehead);
save(filenameStats,'out');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% optional extra stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% count cac in each cluster
% brp is counted by dbscan
% search_radius = 1.5;
% 
% cac = [r.xc(cac_idx) r.yc(cac_idx)];
% brp = [r.xc(brp_idx) r.yc(brp_idx)];
% cac_frame_out = r.frame.*0;
% dist_mat_cac = l2_dist_mat(centers_brp',cac');
% dist_mat_brp = l2_dist_mat(centers_brp',brp');
%%
% cac_count=zeros(length(centers_brp),1);
% brp_count=zeros(length(centers_brp),1);
% for i=1:length(centers_brp)
% %     brp_idx_now = find(dist_mat_brp(:,i)<search_radius^2);
%     cac_idx_now = find(dist_mat_cac(:,i)<search_radius^2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     debugging - show selected molecules based on search radius
% %     imshow(imadjust(count_copy_filter))
% %     hold on
% %     plot(brp(:,1),brp(:,2),'g.')
% %     
% %     axis equal
% %     plot(centers_brp(i,1),centers_brp(i,2),'m+');
% %     plot(brp(brp_idx_now,1),brp(brp_idx_now,2),'b.')
%     cac_count(i) = numel(cac_idx_now);
% %     brp_count(i) = numel(brp_idx_now);
% % keyboard
%     
% 
% end
% %% testing to make sure indexing is correct
% clf
% plot(centers_heatmap(match_idx_filter(:,2),1),centers_heatmap(match_idx_filter(:,2),2),'m+')
% hold on
% axis equal
% plot(centers_brp(match_idx_filter(:,1),1),centers_brp(match_idx_filter(:,1),2),'b+')


% out.cac_count_match = cac_count(match_idx_filter(:,1));
% out.brp_count_match = brp_count(match_idx_filter(:,1));


%% structural averaging

% find Cac within radius around brp cluster centers
storm_mat = StructToMat(r);
storm_mat_filter = storm_mat(brp_idx(brp_idx_filter),:);
heatmap_mat = StructToMat(CreateMolListStruct(hm_x_use,hm_y_use));
%delete non clustered locs

idx_delete = find(cluster_idx_heatmap==0);
heatmap_mat(idx_delete,:)=[];
cluster_idx_heatmap_delete = cluster_idx_heatmap;
cluster_idx_heatmap_delete(idx_delete)=[];

%index heatmap list to corresponding brp cluster
%find matching list by selecting hm cluster idx == matched cluster idx
idx_match_hm = ismember(cluster_idx_heatmap_delete,idx_hm_match_original);
heatmap_mat_match = heatmap_mat(idx_match_hm,:);
cluster_idx_heatmap_match = cluster_idx_heatmap_delete(idx_match_hm);
cluster_idx_heatmap_match_to_brp = cluster_idx_heatmap_match.*0;

for i = 1:length(match_idx)
    idx_change_now = find(cluster_idx_heatmap_match==idx_hm_match_original(i));
    %this sets each heatmap loc to the appropriate frame for Brp alignment
    cluster_idx_heatmap_match_to_brp(idx_change_now) = idx_storm_match_original(i);
end

% for cropping based on clustered data
[heatmap_match_out] = CropROIHR_nmj(heatmap_mat_match,centers_heatmap,...
                                    cluster_idx_heatmap_match,...
                                    cluster_idx_heatmap_match_to_brp);
                                

                                
%set hm cat = 0
heatmap_match_out(:,1)=0;
% %% test
% for i = 1:length(centers_heatmap)
%     idx_now = find(cluster_idx_heatmap_match==i);
%     plot(heatmap_mat_match(idx_now,2),heatmap_mat_match(idx_now,3),'k.')
%     keyboard
% end


%%
filename_out = sprintf('%s-CropROIs.bin',filehead);

cluster_idx_cac = cac_idx.*0;

cac = [r.xc(cac_idx) r.yc(cac_idx)];
dist_mat_cac = l2_dist_mat(centers_brp',cac');

cac_mat = storm_mat(cac_idx,:);

search_radius = 1.5;
%need to append cac list and cac indices corresponding to brp clusters
cac_count = zeros(length(centers_brp),1);
for i = 1:length(centers_brp)
    cac_idx_now = find(dist_mat_cac(:,i)<search_radius^2);
    %set to cluster_idx
    cluster_idx_cac(cac_idx_now)=i;
    cac_count(i) = numel(cac_idx_now);
end
%cat matrices
cluster_idx_combined = cat(1,cluster_idx_brp,cluster_idx_cac);
storm_mat_combined = cat(1,storm_mat_filter,cac_mat);

% CropROIHR_nmj(storm_mat_filter,centers_brp,cluster_idx_brp,filename_out);
% storm_mat_combined = CropROIHR_nmj(storm_mat_combined,centers_brp,cluster_idx_combined,0);
% do this one for replication based cropping

storm_mat_combined = CropROIHR_nmj_replicate(storm_mat_combined,centers_brp,cluster_idx_combined,0);

plot(storm_mat_filter(:,2),storm_mat_filter(:,3),'k.')

mol_list_out = MatToStruct(cat(1,storm_mat_combined,heatmap_match_out));
WriteMolBinNXcYcZc(mol_list_out,filename_out);

%% save output
%
% out.heatmap_count_match = counts_heatmap(match_idx_filter(:,2));
% out.brp_count_match = counts_brp(match_idx_filter(:,1))';
% out.heatmap_centers_match = centers_heatmap(match_idx_filter(:,2),:);
% out.brp_centers_match = centers_brp(match_idx_filter(:,1),:);
% out.heatmap_count_match = out.heatmap_count_match';
% out.cac_count_match = cac_count(match_idx_filter(:,2));
% out.cac_nomatch_idx = find(~ismember((1:length(centers_brp)),match_idx_filter(:,1)));
% out.hm_nomatch_idx = find(~ismember((1:length(centers_heatmap)),match_idx_filter(:,2)));
% out.cac_count_nomatch = cac_count(out.cac_nomatch_idx);
% out.brp_count_nomatch = counts_brp(out.cac_nomatch_idx)';
% out.brp_centers_nomatch = centers_brp(out.cac_nomatch_idx);
% out.hm_centers_nomatch = centers_heatmap(out.hm_nomatch_idx);
% out.hm_counts_nomatch = counts_heatmap(out.hm_nomatch_idx);
% 
% out.match_idx_filter=match_idx_filter;
% out.brp_centers=centers_brp;
% out.heatmap_centers=centers_heatmap;

%% select ROIs for analysis
% - select ROI
% find clusters within ROI
% only select matched clusters to simplify
% analyze_ROIs=1;
% while ~~analyze_ROIs
%     
    clf
    plot(out.heatmap_centers_match(:,1),out.heatmap_centers_match(:,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(out.brp_centers_match(:,1),out.brp_centers_match(:,2),'b+','MarkerSize',5)
    ax=gca;
    
    try
        rect=getrect(ax)
    catch
        analyze_ROIs=0;
%         break
    
    end
    %find points within rectangle
    hm_match = find(out.heatmap_centers_match(:,1)>rect(1)&...
                out.heatmap_centers_match(:,2)>rect(2)&...
                out.heatmap_centers_match(:,1)<(rect(1)+rect(3))&...
                out.heatmap_centers_match(:,2)<(rect(2)+rect(4)));

    brp_match = find(out.brp_centers_match(:,1)>rect(1)&...
                out.brp_centers_match(:,2)>rect(2)&...
                out.brp_centers_match(:,1)<(rect(1)+rect(3))&...
                out.brp_centers_match(:,2)<(rect(2)+rect(4)));
            
            hm_match_2 = match_idx_filter(brp_match,2);
            
    clf
    plot(centers_heatmap(hm_match_2,1),centers_heatmap(hm_match_2,2),'m+','MarkerSize',5)
    hold on
    plot(out.brp_centers_match(brp_match,1),out.brp_centers_match(brp_match,2),'b+','MarkerSize',5)
    
    % find counts and plot
    
    hm_counts_match = counts_heatmap(hm_match_2);
    brp_counts_match = out.brp_count_match(brp_match);
    
    clf
    plot(brp_counts_match,hm_counts_match,'k.')
%    
% %     keyboard
%     
% end

