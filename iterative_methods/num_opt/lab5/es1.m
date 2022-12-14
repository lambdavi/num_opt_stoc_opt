clear all
close all
clc

load("test_functions2.mat")
[xk,fk, gradfk_norm, k, xseq]=newton(x0, f3, gradf3, Hessf3, alpha, kmax, tolgrad);

[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
Z = X.^2 + 4*Y.^2 + 5;
fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X,Y,Z);
hold on
plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
hold off

fig2 = figure();
surf(X, Y, Z, 'EdgeColor', 'none')
hold on
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f1(x0), f1(xseq)], 'r--*')
hold off
