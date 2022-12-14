clear all
close all
clc
load("test_functions2.mat")
load("forcing_terms.mat")
h=sqrt(eps);
c1=1e-4;
rho=0.8;
pcg_maxit=50;
btmax=50;
disp(c1);
[xk, fk, gradfk_norm, k, xseq, btseq]=i_newton_general(x0, f2, gradf2, Hessf2, kmax, tolgrad, c1, rho, btmax, ...
    "", "fw", h, pcg_maxit, fterms_quad);

disp(k);