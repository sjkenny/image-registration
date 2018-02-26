function varname = createUniqueName(baseString)
% Easy creation of unique variable name (evaluated in base workspace)
%
% Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 01/07/2015

% Copyright The MathWorks, Inc. 2015
% Modifications: 
% 6/2/2015: Returns the baseString itself if that is a valid variable name

tmp = evalin('base',['exist(''', baseString,''',''var'')']);
if ~tmp
	varname = baseString;
	return
end
	
n = 0; tmp = 1;
while tmp
	n = n + 1;
	tmp = evalin('base',['exist(''', baseString, num2str(n),''',''var'')']);
end
varname = [baseString, num2str(n)];