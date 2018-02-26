addpath(genpath('C:\STORM'))

close all

R1 = OpenMolList;
R2 = OpenMolList;   %heatmap

X1 = [R1.x, R1.y]';
X2 = [R2.x, R2.y]';

D = l2_dist_mat(X1,X2);

%compute optimal matching
[ids,E_munkres] = munkres(D);

%Choose subset of matches to show
% p  = randperm(N);
% num_shown = 45;
% p = p(1:num_shown);

%shift X2 (for visualization)
X2(2,:) = X2(2,:) + 30;
count = 0;
%Show matching data
figure
hold on;

X3 = X1;
X3(:,ids) = [];
plotp(X1(:,ids), 'b');
plotp(X2, 'r');
plotp(X3,'m');
h = plot([X1(1,ids);X2(1,:)], [X1(2,ids);X2(2,:)], 'k');
set(h, 'LineWidth', 2);
axis('off'); axis('equal');
% title('Matching data');
title(sprintf('%d',count));
set(gcf,'Position',[100 100 900 600]);
file_out = sprintf('s\\outpng_%d.png',count);
saveas(gcf,file_out)
hold off
E_list=E_munkres;

WarpPoints=1;
while WarpPoints==1
    count=count+1;
    x2tform = X2';
    x1tform = X1(:,ids)';
    tform = fitgeotrans(x2tform,x1tform,'projective');
    X2=transformPointsForward(tform,X2');
    X2 = X2';
    D = l2_dist_mat(X1,X2);
    [ids,E_munkres] = munkres(D);
    E_list = cat(1,E_list,E_munkres);
    
    X2(2,:) = X2(2,:) + 30;

    figure
    hold on;

    X3 = X1;
    X3(:,ids) = [];
    plotp(X1(:,ids), 'b');
    plotp(X2, 'r');
    plotp(X3,'m');
    h = plot([X1(1,ids);X2(1,:)], [X1(2,ids);X2(2,:)], 'k');
    set(h, 'LineWidth', 2);
    axis('off'); axis('equal');
    title(sprintf('%d',count));
    set(gcf,'Position',[100 100 900 600]);
    hold off
    keyboard
    file_out = sprintf('s\\outpng_%d.png',count);
    saveas(gcf,file_out)
    close all
end








