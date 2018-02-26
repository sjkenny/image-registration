[FileName,PathName] = uigetfile('*.tif','Select the tif file to convert');
FullFileName=sprintf('%s%s',PathName,FileName);
filehead = FullFileName(1:end-4); %File name without ".bin" or "*.txt"
OutDaxFile = sprintf('%s.dax',filehead);
OutInfFile = sprintf('%s.inf',filehead);

orthophotoRead = imread(FullFileName);
orthophoto=orthophotoRead(:,:,1);
[DimY,DimX]=size(orthophoto);
orthophoto=flipdim(orthophoto,2);

figure, imshow(orthophoto)

%----Output dax-------
fid = fopen(OutDaxFile, 'w');
fwrite(fid, orthophoto', 'uint16','b');
fclose(fid);
fid = fopen(OutInfFile, 'w');
fprintf(fid, 'vstart=1\nhstart=1\nhbin=1\nvbin=1\nnumber of frames = 8000\n');
fprintf(fid, 'vend=%d\n',DimY);
fprintf(fid, 'hend=%d\n',DimX);
fprintf(fid, 'frame size = %d\n',DimX*DimY);
fclose(fid);
%----Output dax-------

