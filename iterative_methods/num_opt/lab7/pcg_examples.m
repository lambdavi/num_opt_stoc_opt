%% DATA INITIALIZATION

clear
clc

n = 10;

A = rand(n);
b = sum(A, 2);

Ax_prod = @(x) A * x;

x_star = ones(n, 1);

tol = 1e-6;
max_iter = 100;

% HELP OF PCG
help pcg

%% PCG: A, b EXAMPLE (PRINTED FLAG)

clc

x_hat = pcg(A, b, tol, max_iter);

disp('PCG Solution:')
x_hat


%% PCG: A, b EXAMPLE (STORED FLAG)

clc

[x_hat, flag] = pcg(A, b, tol, max_iter);


disp(['PCG Solution (flag ', num2str(flag), '):'])
x_hat

%% PCG: Ax, b EXAMPLE (STORED FLAG)

clc

[x_hat, flag] = pcg(Ax_prod, b, tol, max_iter);


disp(['PCG Solution (flag ', num2str(flag), '):'])
x_hat




