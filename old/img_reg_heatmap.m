% total image alignment
% step 1: CP mapping from conventional images
% step 2: CP mapping from rough aligned heatmap+STORM image
% step 3: demon transform
% save tform and displacement field as outputs

% load images
addpath ../common

if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.png',LastFolder);
else
    GetFileName='*.png';
end

[FileName,PathName] = uigetfile(GetFileName,'Select conventional image');
LastFolder=PathName;
GetFileName=sprintf('%s/*.png',LastFolder);

[FileName_s,PathName_s] = uigetfile(GetFileName,'Select STORM image');
LastFolder=PathName;

FullFileName_s=sprintf('%s%s',PathName_s,FileName_s);

FullFileName=sprintf('%s%s',PathName,FileName);
filehead = FullFileName(1:end-4); %File name without ".bin" or "*.txt"
OutDaxFile = sprintf('%s.dax',filehead);
OutInfFile = sprintf('%s.inf',filehead);

filename_tform = sprintf('%s_tform.mat',filehead);


Img = imread(FullFileName);
Img_s = imread(FullFileName_s);
%load in heatmap data

unregistered = ReferenceImage_Crop;
unregistered_adj = imadjust(unregistered);
orthophoto_adj = imadjust(uint16(rgb2gray(Img)));
cpreg=1;
while cpreg==1
    [registered_adj,tform] = ImageRegCP(unregistered_adj,orthophoto_adj);
    imshowpair(registered_adj,orthophoto_adj)
    prompt = 'Check CP mapping (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
end
%%
s = size(orthophoto_adj);
B = imtransform(Probability_Heatmap_Filtered,tform,'XData',...
            [1 s(2)], 'YData', [1 s(1)],...
            'FillValues',0, 'Size', [s(1) s(2)]);
        
        
    idx = find(isnan(B));          

    D=B;
    D(idx)=0;  
    D=imadjust(D);

cpreg=1;
orthophoto_adj_s = imadjust(uint16(rgb2gray(Img_s)));


while cpreg==1
    [registered_adj_s,tform_s] = ImageRegCP(D,Img_s);
    imshowpair(registered_adj_s,orthophoto_adj_s)
    prompt = 'Check CP mapping (press 0 to exit, 1 to retry): ';
    thresh = input(prompt);
    if thresh==0
        cpreg=0;
        break
    end
end

while KeepGoing==1
    prompt = 'Set intensity threshold for heatmap (press 0 to exit): ';

    thresh = input(prompt);
    if thresh==0
        KeepGoing=0;
        break
    end

    thresh_low = thresh/thresh_high;
    if filter_img==1
        img_bw = im2bw(img_subtract,thresh_low);
    else
        img_bw = im2bw(img,thresh_low);
    end
    imshow(img_bw)
end
%% gaussian blur
H=fspecial('gaussian',20,1);
orthophoto_adj_s_blur = imadjust(imfilter(orthophoto_adj_s,H));
imshowpair(registered_adj_s,orthophoto_adj_s_blur)
hold on
%% demon transform
KeepGoing=1;
while KeepGoing==1
    prompt = 'Set field smoothing (press 0 to exit): ';
    FieldSmoothing = input(prompt);
    if ~FieldSmoothing
        KeepGoing=0;
        break
    end   

    [dField,moving_reg] = imregdemons(orthophoto_adj_s_blur,registered_adj_s,[400 200 100],...
        'PyramidLevels',3,'AccumulatedFieldSmoothing',FieldSmoothing);
    movingRegistered = imwarp(orthophoto_adj_s_blur,dField);
    imshowpair(registered_adj_s,movingRegistered)
end

save(filename_tform,'dField','tform');

%% warp heatmap coords
x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];
cat_out = ones(size(x1));
cat_out(length(All_ReMapped_IbIs_Location_Coords_Ib)+1:end)=2;

[bx1,by1] = tformfwd(tform,x1,y1);
[bx2,by2] = tformfwd(tform_s,bx1,by1);
%crop molecule list
img_size = size(Img);
idx_crop = find(bx1>2&bx1<img_size(2)&by1>2&by1<img_size(1));
bx1 = bx1(idx_crop);
by1 = by1(idx_crop);
[cx1,cy1] = WarpPointsDF(bx1,by1,dField);

[r_out] = CreateMolListStruct(cx1,cy1);

OutFile_heatmap=sprintf('%sheatmap_registered.bin',filehead);

WriteMolBinNXcYcZc(r_out,OutFile_heatmap);


