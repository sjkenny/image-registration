% do affine registration as 1st step for point matching

% load images
addpath ../common
addpath old
if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.png',LastFolder);
else
    GetFileName='*.png';
end

[FileName,PathName] = uigetfile(GetFileName,'Select image');
LastFolder=PathName;

FullFileName=sprintf('%s%s',PathName,FileName);
filehead = FullFileName(1:end-4); %File name without ".bin" or "*.txt"
OutDaxFile = sprintf('%s.dax',filehead);
OutInfFile = sprintf('%s.inf',filehead);

filename_tform = sprintf('%s_tform.mat',filehead);


Img = imread(FullFileName);

%load in heatmap data

unregistered = Basal_SynapGCaMP6f_Fluorescence;
unregistered_adj = imadjust(unregistered);
orthophoto_adj = imadjust(uint16(rgb2gray(Img)));
s = size(orthophoto_adj);

cpreg=1;
while cpreg==1
    [registered_adj,tform] = ImageRegCP(unregistered_adj,orthophoto_adj,'similarity');
    imshowpair(registered_adj,orthophoto_adj)
    prompt = 'Check CP mapping (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
end

%warp heatmap coords
%%
% x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
% y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];

%for new data

x1=[Evoked_Coordinates_Ib(:,2); Evoked_Coordinates_Is(:,2);...
    Spont_Coordinates_Ib(:,2); Spont_Coordinates_Is(:,2)];
y1=[Evoked_Coordinates_Ib(:,1); Evoked_Coordinates_Is(:,1);...
    Spont_Coordinates_Ib(:,1); Spont_Coordinates_Is(:,1)];

num_evoked_Ib = length(Evoked_Coordinates_Ib);
num_evoked_Is = length(Evoked_Coordinates_Is);
num_spont_Ib = length(Spont_Coordinates_Ib);
num_spont_Is = length(Spont_Coordinates_Is);

%index each category of coordinates
coords_idx = ones(length(x1),1);
coords_idx(num_evoked_Ib+1:num_evoked_Ib+num_evoked_Is)=2;
coords_idx(num_evoked_Ib+num_evoked_Is+1:num_evoked_Ib+num_evoked_Is+num_spont_Ib)=3;
coords_idx(num_evoked_Ib+num_evoked_Is+num_spont_Ib+1:num_evoked_Ib+num_evoked_Is+num_spont_Ib+num_spont_Is)=4;

cat_out = ones(size(x1));

%cat Ib = 2
% cat_out(length(All_ReMapped_IbIs_Location_Coords_Ib)+1:end)=2;

[bx1,by1] = tformfwd(tform,x1,y1);

%to write to .mat file

%crop molecule list
img_size = size(Img);
idx_crop = find(bx1>2&bx1<img_size(2)&by1>2&by1<img_size(1));
bx1 = bx1(idx_crop);
by1 = by1(idx_crop);
coords_idx_crop = coords_idx(idx_crop);

reg_evoked_Ib = [bx1(coords_idx_crop==1) by1(coords_idx_crop==1)];
reg_evoked_Is = [bx1(coords_idx_crop==2) by1(coords_idx_crop==2)];
reg_spont_Ib = [bx1(coords_idx_crop==3) by1(coords_idx_crop==3)];
reg_spont_Is = [bx1(coords_idx_crop==4) by1(coords_idx_crop==4)];

filename_out=sprintf('%s_heatmap_coords_registered.mat',filehead)
save(filename_out,'reg_evoked_Ib','reg_evoked_Is','reg_spont_Ib','reg_spont_Is','bx1','by1','coords_idx_crop','tform');

cat_out = cat_out(idx_crop);








%%
clf
imshow(Img);
hold on
plot(bx1,by1,'m.')

%%
% 
% [r_out] = CreateMolListStruct(bx1,by1);
% 
% r_out.cat=cat_out;
% 
% OutFile_heatmap=sprintf('%sheatmap_registered_similarity.bin',filehead);
% WriteMolBinNXcYcZc(r_out,OutFile_heatmap);
