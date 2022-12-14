clear all
close all
clc

load("test_functions2.mat")
c1=1e-4;
rho=0.8;
btmax=50;
x0=[-1;-2];
%% FUNCTION 1
[xk,fk, gradfk_norm, k, xseq, btseq]=mod_newton_bcktrck(x0, f1, gradf1, Hessf1, kmax, tolgrad, c1, rho, btmax);

[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
Z = X.^2 + 4*Y.^2 + 5;

fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X,Y,Z);
hold on
plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
hold off
% barplot of btseq
fig2 = figure();
bar(btseq)

fig3 = figure();
surf(X, Y, Z, 'EdgeColor', 'none')
hold on
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f1(x0), f1(xseq)], 'r--*')
hold off

%% FUNCTION 2
[xk,fk, gradfk_norm, k, xseq, btseq]=mod_newton_bcktrck(x0, f2, gradf2, Hessf2, kmax, tolgrad, c1, rho, btmax);

[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
Z = 100*(Y-X.^2).^2 + (1-X).^2;

fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X,Y,Z);
hold on
plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
hold off

% barplot of btseq
fig2 = figure();
bar(btseq)

fig3 = figure();
surf(X, Y, Z, 'EdgeColor', 'none')
hold on
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f2(x0), f2(xseq)], 'r--*')
hold off

%% FUNCTION 3
[xk,fk, gradfk_norm, k, xseq, btseq]=mod_newton_bcktrck(x0, f3, gradf3, Hessf3, kmax, tolgrad, c1, rho, btmax);

[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
Z = (X.^2 +Y -11).^2 + (X+Y.^2 -7).^2;

fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X,Y,Z);
hold on
plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
hold off

% barplot of btseq
fig2 = figure();
bar(btseq)

fig3 = figure();
surf(X, Y, Z, 'EdgeColor', 'none')
hold on
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f3(x0), f3(xseq)], 'r--*')
hold off
