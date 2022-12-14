clear all
close all
clc

load("test_functions2.mat")
h=sqrt(eps);
c1=1e-4;
rho=0.8;
btmax=50;
disp(c1);
[xk, fk, gradfk_norm, k, xseq, btseq]=newton_general(x0, f1, gradf1, Hessf1, kmax, tolgrad, c1, rho, btmax, ...
    "std", "MF", h, 50