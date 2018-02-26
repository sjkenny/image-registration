% 	printf( "PointMatchDemo2 [-g] [-a] [-r] [-f] [-t T_Init] [-e E_Ave] [-i I_Max] fnPntModel fnPntDeform\n" );
% 	printf( "Output of this program is a file with suffix *.match.\n" );
% 	printf( "  -g        -- Do not use relaxation labeling based graph matching.\n" );
% 	printf( "  -a        -- Do not use the affine transformation for the first iteration.\n" );
% 	printf( "  -r        -- Rotation invariant matching.\n" );
% 	printf( "               Avoid using it unless necessary. Imposing rotation invariance\n" );
% 	printf( "               may deteriorate the performance.\n" );
% 	printf( "  -f        -- Force to find as many matches as possible.\n" );
% 	printf( "  -t T_Init -- Parameter related to converting the shape context distance\n" );
% 	printf( "               to a probability measure.\n" );
% 	printf( "               Default: 0.1\n" );
% 	printf( "  -e E_Ave  -- Average number of neighbors of a point.\n" );
% 	printf( "               Default: 5\n" );
% 	printf( "  -i I_Max  -- Maximum number of iterations.\n" );
% 	printf( "               Default: 10\n" );
% 	printf( "  fnPntModel  -- File containing the model point set.\n" );
% 	printf( "  fnPntDeform -- File containing the deformed point set.\n" );
% 	printf( "NOTE: After each iteration of point matching, the model point set is\n" );
% 	printf( "      warped toward the deformed point set. Therefore, the matching procedure\n" );
% 	printf( "      is not symmetric. The point matching results of the following commands may be different.\n" );
% 	printf( "      'PointMatchDemo2 fnPoint1 fnPoint2'\n" );
% 	printf( "      'PointMatchDemo2 fnPoint2 fnPoint1'\n" );


function [match_idx] = pointMatch(points1,points2)
                                    %storm %hm


write_Txt(points1,points2);

exe_name = 'PointMatchDemo_outlier.exe';
%     system(sprintf('%s -a -t %2.2f -f  -i 10 storm_coords.point heatmap_coords.point',exe_name,i))


system(sprintf('%s -g -a -t 0.1 -f -e 5 -i 5 storm_coords.point heatmap_coords.point',exe_name))

% match centers
A=importdata('storm_coords_heatmap_coords.match');
match_idx = A.data;
match_idx=match_idx+1;


clf
plot(points1(:,1),points1(:,2),'m+','MarkerSize',5)
hold on
axis equal
plot(points2(:,1),points2(:,2),'b+','MarkerSize',5)

outlier_thresh=1500;
for j=1:length(match_idx)
    heatmap_idx_now = match_idx(j,2);
    storm_idx_now = match_idx(j,1);
    dist(j) = l2_dist_mat(points2(heatmap_idx_now,:)',points1(storm_idx_now,:)');
    x=[points2(heatmap_idx_now,1) points1(storm_idx_now,1)];
    y=[points2(heatmap_idx_now,2) points1(storm_idx_now,2)];
    if dist(j)>outlier_thresh
        line(x,y,'Color','red')
    else
        line(x,y)
    end
end