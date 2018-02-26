% load image

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

%%
img = imread(FullFileName);
img_gray = (rgb2gray(img));

f1=fspecial('sobel');
f2=fspecial('prewitt');
f3=fspecial('gaussian',5);
f4=fspecial('disk',15);
%%

imf1 = medfilt2(img_gray,[8 8]);
imf2 = medfilt2(imf1,[15 15]);
imf3 = medfilt2(imf1,[20 20]);
imf4 = imbinarize(imf3);

bw=edge(imf4);
imshow(bw)
% imf2 = imfilter(imf1,f2);
% imf3 = imfilter(imf3,f1);
% %%
% img_double = double(imf1);
% img_adjust = img_double./max(img_double(:));
% img_bw = im2bw(img_adjust,0.1);
% 
% imshow(img_bw)



% 
% imf4 = double(imfilter(imf1,f4));
% img_divide = imf1./(0.5.*imf4);
% imshow(imadjust(img_divide))
% imf3 = imfilter(img,f3);
% imshow(imf3)
% 
% imf1 = imfilter(img,f1);
% imf2 = imfilter(img,f2);
 img_bw = im2bw(img_gray,0.1);
