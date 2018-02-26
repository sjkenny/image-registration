%write coords txt files
function write_Txt(storm_coords,heatmap_coords)

outfile_storm=sprintf('storm_coords.point');
outfile_heatmap=sprintf('heatmap_coords.point');

idx = 1:length(heatmap_coords);
idx=(idx-1)';
out_array = [idx heatmap_coords]';


fileID = fopen(outfile_heatmap,'w');
fprintf(fileID,'Point=%d\n',length(heatmap_coords));
fprintf(fileID,'%d %12.8f %12.8f\n',out_array);
fclose(fileID);

idx = 1:length(storm_coords);
idx=(idx-1)';
out_array = [idx storm_coords]';


fileID = fopen(outfile_storm,'w');
fprintf(fileID,'Point=%d\n',length(storm_coords));
fprintf(fileID,'%d %12.8f %12.8f\n',out_array);
fclose(fileID);