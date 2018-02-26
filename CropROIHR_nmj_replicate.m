%this is the version for cropping ROIs of STORM image using coordinates of
%centers from other methods (uses variable 'coords' from HistTxt2D,
%otherwise explicity declare 'coords' as [X;Y] Nx2 matrix before running)
%also crops and aligns .dax file
%for heatmap alignment to STORM, use matched indices to heatmap clusters
%% this one uses set size for ROI and duplicates molecule list


function [storm_mat_out] = CropROIHR_nmj_replicate(storm_mat,centers,cluster_idx,frame_idx)
addpath ../common
offset=10;
roi_size = 5;
ROIHalf=roi_size/2;


% idx_zeros = find(cluster_idx==0);
% cluster_idx(idx_zeros)=[];
% storm_mat(idx_zeros,:)=[];
% 
% cluster_subtract = centers(cluster_idx,:)-offset;
% 
% storm_mat_out = storm_mat;
% %make xy into non-subtracted xc/yc
storm_mat(:,2:3) = storm_mat(:,5:6);

% storm_mat_out(:,5:6)=xy;
% %frame = 15

storm_mat_out = [];
for i=1:length(centers)
    
    cluster_subtract = centers(i,:)-offset;
    fprintf('Cropping %d of %d\n',i,length(centers))
    XMax=centers(i,1)+ROIHalf;
    XMin=centers(i,1)-ROIHalf;
    YMax=centers(i,2)+ROIHalf;
    YMin=centers(i,2)-ROIHalf;
    ROIInd=find(storm_mat(:,2)>XMin&storm_mat(:,2)<XMax&storm_mat(:,3)>YMin&storm_mat(:,3)<YMax);
    data_now = storm_mat(ROIInd,:);
    
    data_now(:,5:6) = data_now(:,5:6)-cluster_subtract;

    data_now(:,15)=i;
    storm_mat_out = cat(1,storm_mat_out,data_now);
    
end
% filename_out='test1.bin';
% WriteMolBinNXcYcZc(MatToStruct(storm_mat_out),filename_out);
