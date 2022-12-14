% Lab 2: Steepest descent
clear
load('test_functions.mat')
f = @(x) sin(0.25*pi*sum(x.^2));
dlgradf = @(x) dlgradient(f(x), x, 'EnableHigherDerivatives', true);
%dlgradf = @(x) dlgradient(f(x), x);

gradf = @(x) dlfeval(dlgradf, x);

[xk, fk, gradfk_norm, k, xseq]=steepest_descent(dlarray(x0), f, gradf, alpha, kmax, tolgrad);