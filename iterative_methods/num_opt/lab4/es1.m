clear
close all

load("test_functions.mat")
alpha0=5;
c1=1e-4;
rho=0.8;
btmax=50;

[xk,fk, gradfk_norm, k, xseq, btseq]=steepest_desc_backtrack(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax);

%[X,Y]=f_meshgrid(xseq, fk');