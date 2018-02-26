%use cpselect to infer transform between images
function [B,mytform] = ImageRegCP(unregistered_adj,orthophoto_adj,method)
s=size(orthophoto_adj);

[xyinput_out, xybase_out] = cpselect(unregistered_adj, orthophoto_adj, 'Wait', true);

mytform=cp2tform(xyinput_out, xybase_out, method);
B = imtransform(unregistered_adj,mytform,'XData',...
                [1 s(2)], 'YData', [1 s(1)],...
                'FillValues',0, 'Size', [s(1) s(2)]);


