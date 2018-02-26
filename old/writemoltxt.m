function f = writemoltxt(mol,filename,append);

if (append)
    fid = fopen(filename,'at');
else
    fid = fopen(filename,'wt');
    N = size(mol.x,1);
    fprintf(fid,'Cas%d	X	Y	Xc	Yc	Height	Area	Width	Phi	Ax	BG	I	Frame	Length	Link	Valid	Z	Zc\n',N);
end

A = [ mol.cat   mol.x  mol.y  mol.x  mol.y  mol.h  mol.area  mol.width  mol.phi  mol.Ax  mol.bg  mol.I  mol.frame  mol.length  mol.link  mol.valid  mol.z  mol.z ];
%   1       2   3   4   5   6       7       8       9   10  11  12
% Cas44178	X	Y	Xc	Yc	Height	Area	Width	Phi	Ax	BG	I	
%   13  14      15      16      17  18
% Frame	Length	Link	Valid	Z	Zc

clear mol;

fprintf(fid,'%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n',A');
fclose(fid);

clear A;
