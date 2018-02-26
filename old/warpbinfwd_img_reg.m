%warpbinfwd for transform defined in image registration app
addpath ..

if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.bin',LastFolder);
else
    GetFileName='*.bin';
end

[FileNameR,PathNameR] = uigetfile(GetFileName,'Select bin file');

RightFile =sprintf('%s%s',PathNameR,FileNameR);
LastFolder=PathNameR;
r= readbinfileNXcYcZc(RightFile);
rx=double(r.xc);
ry=double(r.yc);




GetFileName=sprintf('%s/*.mat',LastFolder);

[FileName,PathName] = uigetfile(GetFileName,'Select warp file');
tformfile =sprintf('%s%s',PathName,FileName);
tform=importdata(tformfile);

[tx,ty] = transformPointsForward(tform,rx,ry);


r.xc=tx;
r.yc=ty;

filehead = RightFile(1:end-4);
outfile = sprintf('%s_Warp_fwd.bin',filehead)
WriteMolBinNXcYcZc(r,outfile);
