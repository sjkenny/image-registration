%save heatmap coords to bin file
%load match file
cat_out = 0;
addpath ../common


[r_out] = CreateMolListStruct(out.fcn_loc_coords(:,1),out.fcn_loc_coords(:,2));


r_out.cat=r_out.cat*0;
r_out.xc = r_out.xc+256;

% OutFile_heatmap=sprintf('%sheatmap_registered_similarity.bin',filehead);
WriteMolBinNXcYcZc(r_out,'heatmap1.bin');