% INITIALIZATION
close all; clear; clc;
load("test_functions2.mat")
disp('** PROBLEM :  ROSENBROCK FUNCTION **');
rho = 0.5; c = 1e-4; kmax = 100; tolgrad = 1e-8;
btmax = 50; n_values = [2];
[X, Y] = meshgrid(linspace(-5.12, 5.12, 500));
Z = 100*(Y -X.^2).^2 + (1-X).^2;
disp('**** NEWTON METHOD WITH BACKTRACKING *****');
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    %x0 = 1.1*ones(n, 1);
    x0 = [-1.2; 1];
tic;
[xk, fk, gradfk_norm, k, xseq, btseq, gmres_it] = newton_bcktrck(x0, f2, gradf2, Hessf2, ...
    kmax, tolgrad, c, rho, btmax,1);

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

    figure();
    % Contour plot with curve levels for each point in xseq
    [C1, ~] = contour(X,Y,Z);
    hold on
    plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
    xlabel("X1");
    ylabel("X2");   
    hold off
end

disp('*** STEEPEST DESCENT WITH BACKTRACKING **');
rho = 0.1;
kmax = 100000; alpha0 = 1; n_values = [2];
for j = 1:length(n_values)
    n = n_values(j);
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);
    % generating starting point
    %x0 = 1.2*ones(n, 1);
tic;
[~, fk, gradfk_norm, k, xseq, btseq] = ...
        steepest_desc_bcktrck(x0, f2, gradf2, alpha0, kmax, ...
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
    figure();
    
    % Contour plot with curve levels for each point in xseq
    [C1, ~] = contour(X,Y,Z);
    hold on
    plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
    xlabel("X1");
    ylabel("X2");   
    hold off

end