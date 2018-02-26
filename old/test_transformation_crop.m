%image registration


%load moving image (STORM)

% moving = rgb2gray(cdata);

%load fixed image (storm)

theta = -0.9;
HorizontalTranslation = 280;
VerticalTranslation = -100;
scale = 1.5;
flip = 1;

tform = [scale*cos(theta) -flip*scale*sin(theta) 0;    ...
         scale*sin(theta) flip*scale*cos(theta) 0;     ...
         HorizontalTranslation VerticalTranslation 1];
T = maketform('affine',tform)



% %create mat of zeros
refsize = size(ref);
moving_size = size(moving);
% storm_conv_resize = padarray(storm_conv,(refsize-storm_conv_size));

%do transformations
%             if ~doflip
%                 movingInputImage = flipud(movingInputImage);
%             end
%            
% RotatedImage = imrotate(storm_conv,theta,'bilinear','crop');
% ScaledImage = imresize(RotatedImage,scale);
% 
% registered = imtranslate(ScaledImage,[HorizontalTranslation,VerticalTranslation]);

registered = imtransform(moving,T,'XData', [1 refsize(2)], 'YData', [1 refsize(1)],'FillValues',0, 'Size', refsize);
% registered = imtransform(storm_conv,T,'FillValues',0, 'Size', refsize);
imshowpair(registered,ref);
% imshowpair(registered,ref);

%identify bounds for cropping
%ymin = first nonzero row
%xmin = first nonzero column
%ymax = last nonzero row
%xmax = last nonzero column
% [row,col] = find(registered);
% xmin = min(col);
% xmax = max(col);
% ymin = min(row);
% ymax = max(row);
% registered_crop = imcrop(registered, [xmin ymin xmax ymax]);
% ref_crop = imcrop(ref, [xmin ymin xmax ymax]);
% 
% imshowpair(registered_crop,ref_crop)


%crop both images
