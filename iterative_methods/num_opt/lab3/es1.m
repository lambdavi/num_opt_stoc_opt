clear all
close all
x_ad = dlarray([30;1]);
%[gradval] = dlfeval(@mygradf,x_ad);
f = @(x) sin(0.25*pi*sum(x.^2));
dlgradf = @(x) dlgradient(f(x), x);
gradf = @(x) dlfeval(dlgradf, x);

