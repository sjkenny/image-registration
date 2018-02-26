function [C] = ImageRegCP(unregistered_adj,orthophoto_adj)

% orthophoto = rgb2gray(Img);
% % orthophoto = imread('647.png');
% % 
% s=size(orthophoto);
% 
% % unregistered = imread('moving.png');
% 
% unregistered = ReferenceImage_Crop;
% unregistered_adj = imadjust(unregistered);
% orthophoto_adj = imadjust(orthophoto);

[xyinput_out, xybase_out] = cpselect(unregistered_adj, orthophoto_adj, 'Wait', true);

mytform=cp2tform(xyinput_out, xybase_out, 'projective');
B = imtransform(unregistered_adj,mytform,'XData',...
                [1 s(2)], 'YData', [1 s(1)],...
                'FillValues',0, 'Size', [s(1) s(2)]);
imshowpair(B,orthophoto)

C = imtransform(Probability_Heatmap_Filtered,mytform,'XData',...
                [1 s(2)], 'YData', [1 s(1)],...
                'FillValues',0, 'Size', [s(1) s(2)]);
            
idx = find(isnan(C));          

D=C;
D(idx)=0;
E=uint16(D.*10000);
%write dax


            
            
% 
% 
% 
% 
% 
% bx=xybase_out(:,1);
% by=xybase_out(:,2);
% 
% ix=xyinput_out(:,1);
% iy=xyinput_out(:,2);
% 
% [tx,ty] = tforminv(mytform,bx,by);
% plot(ix,iy,'k.',tx,ty,'m.');
% 
% save ('test1.mat')
% 
% [FileName,PathName] = uigetfile('*.bin','Select the bin file to warp');
% 
% FullFileName=sprintf('%s%s',PathName,FileName);
% filehead = FullFileName(1:end-4); %File name without ".bin" or "*.txt"
% outfile = sprintf('%s_warp.bin',filehead);
% 
% fprintf(1, 'Reading...');
% r = readbinfileNXcYcZc(FullFileName);
% fprintf(1, 'done!\n');
% 
% bx = double(r.xc);
% by = double(r.yc);
% 
% [tx,ty] = tformfwd(mytform,bx,by);
% 
% r.xc=tx;
% r.yc=ty;
% 
% fprintf(1, 'Writing...');
% WriteMolBinNXcYcZc(r,outfile);
% fprintf(1, 'done!\n');
