%this is the version for cropping ROIs of STORM image using coordinates of
%centers from other methods (uses variable 'coords' from HistTxt2D,
%otherwise explicity declare 'coords' as [X;Y] Nx2 matrix before running)
%also crops and aligns .dax file
%for heatmap alignment to STORM, use matched indices to heatmap clusters


function [storm_mat_out] = CropROIHR_nmj(storm_mat,centers,cluster_idx,frame_idx)
addpath ../common
offset=10;


idx_zeros = find(cluster_idx==0);
cluster_idx(idx_zeros)=[];
storm_mat(idx_zeros,:)=[];

cluster_subtract = centers(cluster_idx,:)-offset;

storm_mat_out = storm_mat;
%make xy into non-subtracted xc/yc
storm_mat_out(:,2:3) = storm_mat_out(:,5:6);
xy = storm_mat(:,5:6)-cluster_subtract;
storm_mat_out(:,5:6)=xy;
%frame = 15
if ~frame_idx
    frame_idx=cluster_idx;
end
storm_mat_out(:,15)=frame_idx;

% filename_out='test1.bin';
% WriteMolBinNXcYcZc(MatToStruct(storm_mat_out),filename_out);

