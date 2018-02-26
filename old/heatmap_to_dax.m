%heatmap to dax file for cluster localization

%insert name here
%preprocess

phf = Probability_Heatmap_Filtered;
phf(isnan(phf))=0;
% phf_adj = imadjust(phf);

phf_int = cast(phf.*10000,'uint16');


filename = 'test1';
info = write_info_template(Probability_Heatmap_Filtered,filename);
WriteDAXFiles(phf_int,info);

