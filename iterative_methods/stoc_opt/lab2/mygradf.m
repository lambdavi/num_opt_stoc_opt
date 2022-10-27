function [dydx] = mygradf(x)
%MYGRADF Compute gradient using DL Toolbox and Backwards AD
%   Input: vector X
%   Output: gradient of Y in the points X
y=sin(0.25*pi*sum(x.^2));
dydx = dlgradient(y,x);
end