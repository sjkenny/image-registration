% linear regression
% need column vectors

function [b1,b2,r1,r2] = linreg(x,y)

a = [(ones(size(y))) x];

b1=x\y;
b2=a\y;

yCalc1 = b1.*x;
yCalc2 = a*b2;

r1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
r2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
