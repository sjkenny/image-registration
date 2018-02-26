% total image alignment
% step 1: filter input images
% step 2: do control point based initial transform
% step 3: do demon transform
% save tform and displacement field as outputs

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
%filtering stuff
H = fspecial('disk',30);

imshowpair(registered_filter,orthophoto_filter)
imshowpair(reg_BG,ortho_BG)

registered_filter_adj = imadjust(registered_filter);
orthophoto_filter_adj = imadjust(orthophoto_filter);
imshowpair(registered_filter_adj,orthophoto_filter_adj)

%% do minimum thresholding
thresh=50;
registered_thresh = registered_filter_adj;
orthophoto_thresh = orthophoto_filter_adj;

imshowpair(registered_thresh,orthophoto_thresh)
KeepGoing=1;
%do local intensity normalization
H = fspecial('disk',30);
ortho_2 = imfilter(orthophoto_thresh,H,'replicate');
reg_2 = imfilter(registered_thresh,H,'replicate');
while KeepGoing==1
    prompt = 'Set intensity threshold (press 0 to exit): ';

    thresh = input(prompt);
    if thresh==0
        KeepGoing=0;
        break
    end
    ortho_BG = imfilter(orthophoto_adj,H,'replicate');
    reg_BG = imfilter(registered_adj,H,'replicate');
    registered_filter = registered_adj-thresh.*reg_BG;
    orthophoto_filter = orthophoto_adj-thresh.*ortho_BG;
    registered_filter_adj = imadjust(registered_filter);
    orthophoto_filter_adj = imadjust(orthophoto_filter);
    
    
%     
%     registered_adj_filter = registered_adj;
%     orthophoto_adj_filter = orthophoto_adj;
%     
%     registered_adj_filter(find(reg_2<thresh))=0;
%     orthophoto_adj_filter(find(ortho_2<thresh))=0;
%     
%     registered_adj_bw = im2bw(registered_adj,thresh);
%     orthophoto_adj_bw = im2bw(orthophoto_adj,thresh);

    imshowpair(registered_filter_adj,orthophoto_filter_adj)
    hold on
end


registered_adj_filter(find(reg_2<thresh))=0;
orthophoto_adj_filter(find(ortho_2<thresh))=0;
imshowpair(registered_adj_filter,orthophoto_adj_filter)
% ortho_3 = 0.5.*ortho_2+1;
% orthophoto_divide = bsxfun(@rdivide,orthophoto_thresh,ortho_3);
% 
% reg_3 = 0.5.*reg_2+1;
% registered_divide = bsxfun(@rdivide,registered_thresh,reg_3);
% 
% imshowpair(imadjust(registered_divide),imadjust(orthophoto_divide))

%%
%with filtering
% [dField,moving_reg] = imregdemons(registered_thresh,orthophoto_thresh,[400 200 100],...
%     'PyramidLevels',3,'AccumulatedFieldSmoothing',1.5);
% movingRegistered = imwarp(registered_thresh,dField);
% imshowpair(movingRegistered,orthophoto_thresh)
%without filtering
KeepGoing=1;
while KeepGoing==1
    prompt = 'Set field smoothing (press 0 to exit): ';

    FieldSmoothing = input(prompt);
    if ~FieldSmoothing
        KeepGoing=0;
        break
    end
    
%warp storm to heatmap - dfield will then warp heatmap to storm
    [dField,moving_reg] = imregdemons(orthophoto_filter_adj,registered_filter_adj,[400 200 100],...
        'PyramidLevels',3,'AccumulatedFieldSmoothing',FieldSmoothing);
    movingRegistered = imwarp(orthophoto_adj,dField);
    imshowpair(registered_adj,movingRegistered)
    
end
    
%without filtering

% [dField,moving_reg] = imregdemons(registered_divide,orthophoto_divide,[800 400 200 100],...
%     'PyramidLevels',4,'AccumulatedFieldSmoothing',1);
% movingRegistered = imwarp(registered_divide,dField);  
% imshowpair(movingRegistered,orthophoto_filter_adj)
%%

s = size(orthophoto_adj);
%warp heatmap
B = imtransform(Probability_Heatmap_Filtered,tform,'XData',...
                [1 s(2)], 'YData', [1 s(1)],...
                'FillValues',0, 'Size', [s(1) s(2)]);
            
idx = find(isnan(B));          

D=B;
D(idx)=0;            
C = imwarp(D,dField);
C=imadjust(C);
imshowpair(registered_adj,C)

save(filename_tform,'dField','tform');
%%
%warp heatmap coords
x1=[All_ReMapped_IbIs_Location_Coords_Ib(:,2); All_ReMapped_IbIs_Location_Coords_Is(:,2)];
y1=[All_ReMapped_IbIs_Location_Coords_Ib(:,1); All_ReMapped_IbIs_Location_Coords_Is(:,1)];
cat_out = ones(size(x1));
cat_out(length(All_ReMapped_IbIs_Location_Coords_Ib)+1:end)=2;


[bx1,by1] = tformfwd(tform,x1,y1);


%crop molecule list
img_size = size(Img);
idx_crop = find(bx1>2&bx1<img_size(2)&by1>2&by1<img_size(1));
bx1 = bx1(idx_crop);
by1 = by1(idx_crop);
[cx1,cy1] = WarpPointsDF(bx1,by1,dField);


% 
% [r,filename] = OpenMolList;
% xStorm = r.xc;
% yStorm = r.yc;
% [dx1,dy1] = WarpPointsDF(xStorm,yStorm,dField);

%write bin file

% 
% plot(x1,y1,'k.')
% hold on
% plot(dx1,dy1,'m.')
% hold off

% imshow(dd)
% hold on
% plot(bx1,by1,'m.')
% hold off
%%
imshowpair(registered_filter,movingRegistered)
hold on
plot(cx1,cy1,'g.')

hold off
%%
% %C is green
% imshow(orthophoto_adj)
% hold on
% plot(dx1,dy1,'m.')
% plot(bx1,by1,'m.')

[r_out] = CreateMolListStruct(cx1,cy1);
% 
% r_out_storm = r;
% r_out_storm.xc = dx1;
% r_out_storm.yc = dy1;


% OutFile=sprintf('%s_registered.bin',filehead);
OutFile_heatmap=sprintf('%sheatmap_registered.bin',filehead);
% WriteMolBinNXcYcZc(r_out_storm,OutFile);
WriteMolBinNXcYcZc(r_out,OutFile_heatmap);










