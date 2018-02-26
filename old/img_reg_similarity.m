% do affine registration as 1st step for point matching

% load images
addpath ../common

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

unregistered = ReferenceImage_Crop;
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
x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];
cat_out = ones(size(x1));

%cat Ib = 2
cat_out(length(All_ReMapped_IbIs_Location_Coords_Ib)+1:end)=2;

[bx1,by1] = tformfwd(tform,x1,y1);
%crop molecule list
img_size = size(Img);
idx_crop = find(bx1>2&bx1<img_size(2)&by1>2&by1<img_size(1));
bx1 = bx1(idx_crop);
by1 = by1(idx_crop);

cat_out = cat_out(idx_crop);
%%
clf
imshow(Img);
hold on
plot(bx1,by1,'m.')

%%
[r_out] = CreateMolListStruct(bx1,by1);

r_out.cat=cat_out;

OutFile_heatmap=sprintf('%sheatmap_registered_similarity.bin',filehead);
WriteMolBinNXcYcZc(r_out,OutFile_heatmap);
