function [img_copy center_list idx_list] = flood_fill_indices(img_bw)

%input = binary or grayscale image
%variable input = minimum intensity threshold for connectivity
%output = indexed flood fill matrix
%         center of mass of each area
%         linear indices of each area
%finish thresh_low later

if nargin > 1
    thresh_low = varargin{1};
else
    thresh_low = 0;
end


% flood fill with indices
% 
% test file
% load('flood_fill_test.mat')
%
%drag in flood_fill_test.mat
%target_color = 1
%replace_idx = 0
%use im_bw to check indices, use im_copy for output
imgsize=size(img_bw);
% vr = VideoWriter('floodfill_ind.mp4','MPEG-4');
% open(vr);

% img_copy = double(img_bw);
img_copy = zeros(imgsize);
% imshow(img_copy);
FloodFill = 0;
ind = 1;
q = [];
%idx for each structure
replace_idx = 1;
while ~FloodFill
    %find white pixel
    [i,j] = find(img_bw>0,1,'first');
    if isempty(i)
        FloodFill = 1;
        break
    end
    %remove from list of pixels to check, add node to queue
    q_add = [i,j];
    q = [q; q_add];
%     im_copy(i,j)=1;
%     img_bw(i,j)=0;
    %check all adjacent pixels, add to queue
    while ~isempty(q)
        %move right/left to bounds
        x = q(1,2);
        y = q(1,1);
        
        edge_left = 0;
        edge_right = 0;
        while ~edge_left
            %check for nonzero index
            if x>1
                if img_bw(y,x-1)==0
                    edge_left = x;
                end
            else
                edge_left = x;
            end
            %step left
            x = x-1;
        end
        x = q(1,2);
        while ~edge_right
            if x<imgsize(2)
                if img_bw(y,x+1)==0
                    edge_right = x;
                end
            else
                edge_right = imgsize(2);
            end
            %step right
            x=x+1;
        end
        %fill row, check above/below each node
        img_bw(y,edge_left:edge_right) = 0;
        img_copy(y,edge_left:edge_right) = replace_idx;
        
%         imshow(img_bw)
        for m = edge_left:edge_right
            if y+1<imgsize(1)
                if img_bw(y+1,m)>0
                    q_add = [y+1,m];
                    q = [q;q_add];
                end
            end
            if y-1>0
                if img_bw(y-1,m)>0
                    q_add = [y-1,m];
                    q = [q;q_add];
                end
            end
        end
        %delete node from queue
        q(1,:)=[];
    end

    %find CoM
    
    linear_idx=find(img_copy==replace_idx);
    I=img_bw(linear_idx);
    [row,col]=ind2sub(imgsize,linear_idx);
    idx_list{replace_idx} = linear_idx;
    center_list(replace_idx,1)=sum(col.*linear_idx)/sum(linear_idx);
    center_list(replace_idx,2)=sum(row.*linear_idx)/sum(linear_idx);
%     
    replace_idx = replace_idx+1;
    
        
%         
%     if img_bw(i+1,j+1)==1
%         q_add = [i,j];
%         q = [q;q_add];
        

end

        
    