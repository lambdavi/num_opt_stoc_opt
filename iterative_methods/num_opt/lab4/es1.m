clear all
close all
clc
load("test_functions.mat")
alpha0=5;
c1=1e-4;
rho=0.8;
btmax=50;
disp('STEEPST DESCENT: START')
[xk, fk, gradfk_norm, k, xseq, btseq] = steepest_desc_backtrack(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax);
disp('STEEPST DESCENT: FINISHED')
disp('STEEPST DESCENT: RESULTS')
disp(['xk: ', mat2str(xk), '(actual minimun: [0; 0]);'])
disp(['f(xk): ', num2str(fk), '(actual min. value: 0);'])
disp(['N iterations: ', num2str(k), '/', num2str(kmax)])

%% PLOTS

% Creation of the meshgrid for the counter-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% compute the values of f for each point of the mesh
Z = X.^2 + Y.^2;

% Plots
f1 = figure();
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
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f(x0), f(xseq)], 'r--*')
hold off