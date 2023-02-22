% Script for minimizing the function in problem 25 (extended rosenbrock)
% INITIALIZATION
close all; clear; clc;
disp('** PROBLEM 25: EXTENDED ROSENBROCK FUNCTION **');
rho = 0.5; c = 1e-4; kmax = 100; tolgrad = 1e-8;
btmax = 50; n_values = [1e3, 1e4];
alpha0=1;
% Function handles
f = @(x) problem_25_function(x); % value of the function
gradf = @(x) problem_25_grad(x); % gradient vector
Hessf = @(x) problem_25_hess(x); % hessian matrix

disp('**** NEWTON METHOD WITH BACKTRACKING *****');
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    x0 = zeros(n, 1);
    for i = 1:n
        if mod(i,2) == 1
            x0(i) = -1.2;
        else
            x0(i) = 1.0;
        end 
    end
tic;
[xk, fk, gradfk_norm, k, xseq, btseq, gmres_it, ratio] = newton_bcktrck(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c, rho, btmax,0);

elapsed_time = toc;
disp('************** RESULTS ****************');
disp(['f(xk): ', num2str(fk(end))]);
disp(['gradfk_norm: ', num2str(gradfk_norm(end))]);
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';']);
disp(['ratio: ', num2str(median(ratio))]);
disp(['Elapsed time: ', num2str(elapsed_time, '%.3f') ' sec']);
disp('***************************************');
% line plot of function value and gradient norm
figure();

yyaxis left;
    semilogy(fk, 'LineWidth', 2);
    ylabel('Value of the function');
    yyaxis right;
    semilogy(gradfk_norm, '--', 'LineWidth', 2);
    ylabel('Norm of the gradient');
    title({'[Problem 25]' 'Newton method, n=', num2str(n, '%.0e')});
    legend('Fk trend', 'GradFk trend', 'Location', 'southwest');
    xlabel('Number of iteration (k)');
    disp('---');
    % histogram plot of btseq values
    figure();
    histogram(btseq);
    title({'[Problem 25]' 'Histogram backtracking iterations'...
        'Newton method n=', num2str(n, '%.0e')});
    xlabel('Number of backtracking iterations');
    ylabel('Newton method iterations');
    xticks(0:3);
    % bar plot of gmres number of iterations
    figure();
    bar(1:k, gmres_it);
    ylim([0, 5]);
    title({'[Problem 25]' 'PCG iterations'});
    xlabel('Iterations of the Newton method');
    ylabel('Number of PCG iterations')
end



disp('*** STEEPEST DESCENT WITH BACKTRACKING **');
kmax = 100000; alpha0 = 1; n_values = [1e3, 1e4];
rho=0.2;
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    x0 = zeros(n, 1);
    for i = 1:n
        if mod(i,2) == 1
            x0(i) = -1.2;
        else
            x0(i) = 1.0;

        end 
    end
tic;
[~, fk, gradfk_norm, k, ~, btseq, ratio] = ...
        steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, ...
            tolgrad, c, rho, btmax);
    elapsed_time = toc;
    disp('************** RESULTS ****************');
    disp(['f(xk): ', num2str(fk(end)), ' (actual min. value: 0);']);
    disp(['gradfk_norm: ', num2str(gradfk_norm(end))]);
    disp(['ratio: ', num2str(median(ratio))]);
    disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';']);
    disp(['Elapsed time: ' num2str(elapsed_time, '%.3f') ' sec']);
    disp('***************************************');
    % line plot of function value and gradient norm
    figure();
    yyaxis left;
    semilogy(fk, 'LineWidth', 2);
    ylabel('Value of the function');
    yyaxis right;
    semilogy(gradfk_norm, '--', 'LineWidth', 2);
    ylabel('Norm of the gradient');
    title({'[Problem 25]' 'Gradient method, n=', num2str(n, '%.0e')});
    legend('Fk trend', 'GradFk trend', 'Location', 'northeast');
    xlabel('Number of iteration (k)');
    % histogram plot of btseq values
    figure();
    histogram(btseq);
    title({'[Problem 25]' 'Histogram backtracking iterations,' ...
        'gradient method n=', num2str(n, '%.0e')});
    xlabel('Number of backtracking iterations');
    ylabel('Gradient method iterations');
    xticks(0:10)
end



function Fx = problem_25_function(X)
% Function for computing the value F(x) where F is the extended rosenbrock
%
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% Fx: the value of the extended rosenbrock function at X
f_k_odd = @(x, k) 100 * (x(k).^2 - x(k+1)) .^ 2;
f_k_even = @(x, k) (x(k-1) - 1) .^ 2;
n = length(X);
Fx = 0;
for i=1:n
    if mod(i, 2) == 1
        Fx = Fx + f_k_odd(X, i);
    else
        Fx = Fx + f_k_even(X, i);
    end 
end
Fx = Fx ./ 2;
end

function gradFx = problem_25_grad(X)
% Function for computing the gradient vector in a specified point of the
% extended rosenbrock
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% gradFx: n by 1 vector which represent the gradient at X
grad_odd = @(x, k) 200 * x(k) .^ 3 - 200 .* x(k) .* x(k+1) + x(k) - 1;
grad_even = @(x, k) 100 * x(k) - 100 * x(k-1) .^ 2;
n = length(X);
gradFx = zeros(n, 1);
for i = 1:n
    if mod(i, 2) == 1
        gradFx(i) = grad_odd(X, i);
    else
        gradFx(i) = grad_even(X, i);
    end 
end
end

function HessFx = problem_25_hess(X)
% Function for computing the hessian matrix of extended rosenbrock at a
% given point X 
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% HessFx: a n-by-n sparse tridiagonal matrix representing the hessian
% matrix at point X
    n = length(X);
    max_nz = 2*n - 1;
    row = zeros(max_nz, 1);
    column = zeros(max_nz, 1);
    values = zeros(max_nz, 1);
    nz = 0;
    for i=1:(n-1)
        if mod(i,2) == 1
           nz = nz+1;
           row(nz) = i;
           column(nz) = i;
           values(nz) = 600 * X(i).^2 - 200 * X(i+1) + 1;
           tmp = -200 * X(i);
           nz = nz + 1;
           row(nz) = i;
           column(nz) = i+1;
           values(nz) = tmp;
           nz = nz + 1;
           row(nz) = i+1;
           column(nz) = i;
           values(nz) = tmp;
        else
            nz = nz + 1;
            row(nz) = i;
            column(nz) = i;
            values(nz) = 100;
        end 
    end

    nz = nz + 1;
    row(nz) = n;
    column(nz) = n;
    values(nz) = 100;
    HessFx = sparse(row, column, values, n, n);
end
