addpath ../common
addpath(genpath('image_registration'))
close all

R1 = OpenMolList;
R2 = OpenMolList;

X1 = [R1.x; R1.y];
X2 = [R2.x; R2.y];

R_mat = @(theta) ([cos(theta), -sin(theta); sin(theta), cos(theta)]);

%shift by mean
X1 = 10*bsxfun(@minus,X1,mean(X1,2));
X2 = 10*bsxfun(@minus,X2,mean(X2,2));

% raw estimation of rotation 
X2 = R_mat(-pi/3)*X2;


%Make # points match by subsampling
N = size(X2,2);
p = randperm(size(X1,2));
X1 = X1(:,p(1:N));

%Show raw data
plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

figure
hold on;
plotp(X1, 'b');
plotp(X2, 'r');
axis('off'); axis('equal');
title('Raw points, mean shifted and rotated');

%Compuet pairwise distances efficiently (Although realistically 100 points
%is extremely small)
D = l2_dist_mat(X1,X2);

%compute optimal matching
[ids,E_munkres] = munkres(D);

%Choose subset of matches to show
p  = randperm(N);
num_shown = 45;
p = p(1:num_shown);

%shift X2 (for visualization)
X2(1,:) = X2(1,:);

%Show matching data
figure
hold on;
% h = plot([X1(1,ids(p));X2(1,p)], [X1(2,ids(p));X2(2,p)], 'k');
% set(h, 'LineWidth', 2);
plotp(X1, 'b');
plotp(X2, 'r');
axis('off'); axis('equal');
title('Matching data');












