%plot matches
%load in output
%optional - load in match_idx_crop (saved in quantify_matches)


addpath ../common 

if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.png',LastFolder);
else
    GetFileName='*.png';
end


[FileNameL,PathNameL] = uigetfile(GetFileName,'Select Image 1');

GetFileName=sprintf('%s/*.ini',PathNameL);
[FileNameInf,PathNameInf] = uigetfile(GetFileName,'Select inf file');

ini_name = sprintf('%s%s',PathNameInf,FileNameInf);


a = ReadInfoFile(ini_name);
zoom = a.Zoom;
xStart = a.DispX;
yStart = a.DispY;


%%
fcn_centers = out.fcn_centers;

fcn_centers(:,1) = fcn_centers(:,1)+256;

%coords start halfway through pixel
shift = [xStart-0.5 yStart-0.5];


storm_centers = out.storm_centers;
% fcn_centers = fcn_centers.*zoom;

storm_centers = bsxfun(@minus,storm_centers,shift);
fcn_centers = bsxfun(@minus,fcn_centers,shift);

storm_centers = storm_centers.*2;
fcn_centers = fcn_centers.*2;
match_idx = out.match_idx;


LeftFile = sprintf('%s%s',PathNameL,FileNameL);
LeftImg=imread(LeftFile);
imshow(LeftImg)
hold on
%%
fcn_centers = out.fcn_centers;
storm_centers = out.storm_centers;


fcn_idx_match = match_idx_crop(:,2);
storm_idx_match = match_idx_crop(:,1);
plot(fcn_centers(fcn_idx_match,1),fcn_centers(fcn_idx_match,2),'wo','MarkerSize',5)
% hold on
% axis equal
plot(storm_centers(storm_idx_match,1),storm_centers(storm_idx_match,2),'wo','MarkerSize',5)



for j=1:length(match_idx_crop)
    %referenced back to originial uncropped list
    heatmap_idx_now = match_idx_crop(j,2);
    storm_idx_now = match_idx_crop(j,1);

    x=[fcn_centers(heatmap_idx_now,1) storm_centers(storm_idx_now,1)];
    y=[fcn_centers(heatmap_idx_now,2) storm_centers(storm_idx_now,2)];

    line(x,y)

end

filename_out = sprintf('%s-match.png',LeftFile(1:end-4));
print(gcf,filename_out,'-dpng','-r600');
