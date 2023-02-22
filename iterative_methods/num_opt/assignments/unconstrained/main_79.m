% Script for minimizing the function in problem 25 (extended rosenbrock)
% INITIALIZATION
close all; clear; clc;
disp('** PROBLEM 79: EXTENDED ROSENBROCK FUNCTION **');
rho = 0.5; c = 1e-4; kmax = 10000; tolgrad = 1e-8;
btmax = 50; n_values = [1e3, 1e4];
alpha0=1;
% Function handles

f = @(x) 0.5*sum(small_fx_79(x).^2);
gradf = @(x) problem_79_grad(x);
Hessf = @(x) problem_79_hess(x);


disp('**** NEWTON METHOD WITH BACKTRACKING *****');
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    x0 = -1*ones(n, 1);
    
tic;
[xk, fk, gradfk_norm, k, xseq, btseq, gmres_it, ratio] = newton_bcktrck(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c, rho, btmax,0);

elapsed_time = toc;
disp('************** RESULTS ****************');
disp(['f(xk): ', num2str(fk(end))]);
disp(['gradfk_norm: ', num2str(gradfk_norm(end))]);
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';']);
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
    title({'[Problem 25]' 'GMRES iterations'});
    xlabel('Iterations of the Newton method');
    ylabel('Number of GMRES iterations')
end


disp('*** STEEPEST DESCENT WITH BACKTRACKING **');
kmax = 100000; alpha0 = 1; n_values = [1e3, 1e4];
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    x0 = -1*ones(n, 1);
tic;
[~, fk, gradfk_norm, k, ~, btseq] = ...
        steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, ...
            tolgrad, c, rho, btmax);
    elapsed_time = toc;
    disp('************** RESULTS ****************');
    disp(['f(xk): ', num2str(fk(end)), ' (actual min. value: 0);']);
    disp(['gradfk_norm: ', num2str(gradfk_norm(end))]);
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



%{
disp('*** INEXACT NEWTON WITH BACKTRACKING **');
kmax = 1000; alpha0 = 1;  h=1e-8; n_values=[1e3, 1e4, 1e5];
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
    load("forcing_terms.mat")
    pcg_maxit=50; forcing_terms = fterms_quad;
    [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    i_newton_general(x0, f, gradf, Hessf, kmax, ...
    tolgrad, c, rho, btmax, "", "MF", h, pcg_maxit, forcing_terms);
    elapsed_time = toc;
    disp('************** RESULTS ****************');
    disp(['f(xk): ', num2str(fk(end)), ' (actual min. value: 0);']);
    disp(['gradfk_norm: ', num2str(gradfk_norm(end))]);
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


%}

function small_fx_array = small_fx_79(x)
    n = length(x);
    small_fx_array = zeros(n,1);
    for k = 1:n
        small_fx_array(k) = (3*x(k) - 0.1*x(k).^2 + 1);
        if(k==1)
            small_fx_array(k) = small_fx_array(k) - 2*x(k+1);
        elseif (k==n)
            small_fx_array(k) = small_fx_array(k) - x(k-1);
        else
            small_fx_array(k) = small_fx_array(k) - x(k-1) - 2*x(k+1);
        end
    end
end

function gradFx = problem_79_grad(x)
% Function for computing the gradient vector at a specified point of the
% function reported in problem 76 in [1]
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% gradFx: n by 1 vector which represent the gradient at X
    small_fx_array = small_fx_79(x);
    n = length(x);
    gradFx = zeros(n,1);
    gradFx(1) = (3 - 0.2*x(1))*small_fx_array(1) - small_fx_array(2);
    gradFx(n) = (3 - 0.2*x(n))*small_fx_array(n) - 2*small_fx_array(n-1);
    for i = 2:(n-1)
        gradFx(i) = -2*small_fx_array(i-1) + (3 - 0.2*x(i))*small_fx_array(i) - small_fx_array(i+1);
    end 
end

function HessFx = problem_79_hess(x)
% Function for computing the gradient vector at a specified point of the
% function reported in problem 76 in [1]
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% gradFx: n by 1 vector which represent the gradient at X
    small_fx_array = small_fx_79(x);
    n = length(x);
    
    %Initialize vectors for sparse matrix:
    rows=zeros(5*n-6,1);
    cols=zeros(5*n-6,1);
    vals=zeros(5*n-6,1);
    count=3;
    % Compute the diagonal first and last element
    % HessFx(1,1) = ((3 - 0.2*x(1)).^2)*(-0.2*small_fx_array(1)) + 1;
    rows(1) = 1; cols(1) = 1; vals(1) = ((3 - 0.2*x(1)).^2)*(-0.2*small_fx_array(1)) + 1;
    % HessFx(n,n) = 4 + ((3 - 0.2*x(n)).^2)*(-0.2*small_fx_array(n));
    rows(2) = n; cols(2) = n; vals(2) = 4 + ((3 - 0.2*x(n)).^2)*(-0.2*small_fx_array(n));
    for i = 1:(n-1)
        if(i~=1)
            % Compute the diagonal
            % HessFx(i,i) = 5 + ((3 - 0.2*x(i)).^2)*(-0.2*small_fx_array(i));
            rows(count) = i; cols(count) = i; vals(count) = 5 + ((3 - 0.2*x(i)).^2)*(-0.2*small_fx_array(i));
            count = count+1;
        end
        % Compute the "second" diagonal
        % HessFx(i,i+1) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        rows(count) = i; cols(count) = i+1; vals(count) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        count = count+1;
        % HessFx(i+1,i) = HessFx(i,i+1);
        rows(count) = i+1; cols(count) = i; vals(count) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        count = count+1;

        % Compute the "third" diagonal
        if(i~=n-1)
            % HessFx(i,i+2) = 2;
            rows(count) = i; cols(count) = i+2; vals(count) = 2;
            count = count+1;
            % HessFx(i+2,i) = 2;
            rows(count) = i+2; cols(count) = i; vals(count) = 2;
            count = count+1;
        end
    end 
    HessFx = sparse(rows,cols,vals,n,n);
end