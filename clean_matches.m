function [idx_out] = clean_matches(points1,points2)
%select outliers in points2 (STORM clusters)

outlier_idx=[];
select_rect=1;
idx_out=1:length(points2);
while select_rect==1
    clf
    plot(points1(idx_out,1),points1(idx_out,2),'m+','MarkerSize',5)
    hold on
    axis equal
    plot(points2(idx_out,1),points2(idx_out,2),'b+','MarkerSize',5)
    
    for j=1:length(idx_out)     
        
        x=[points1(idx_out(j),1) points2(idx_out(j),1)];
        y=[points1(idx_out(j),2) points2(idx_out(j),2)];
        line(x,y)
        
    end
    ax=gca;    
    try
        rect=getrect(ax);
    catch
        select_rect=0;
        break
    end
       
    points2_outliers_temp = find(points2(:,1)>rect(1)&...
                                points2(:,2)>rect(2)&...
                                points2(:,1)<(rect(1)+rect(3))&...
                                points2(:,2)<(rect(2)+rect(4)));

    idx_out = 1:length(points2);                         
    outlier_idx = cat(1,outlier_idx,points2_outliers_temp);
    outlier_idx = unique(outlier_idx);
    idx_out(outlier_idx)=[];
    
end

