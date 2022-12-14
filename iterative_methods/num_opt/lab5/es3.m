clear all
close all
clc

load("test_functions2.mat")
c1=1e-4;
rho=0.8;
btmax=50;

f_a = @(x)(-20*exp(-0.2*sqrt(0.5*(x(1,:)^2+x(2,:)^2) - exp(0.5*(cos(2*pi*x(1,:))+cos(2*pi*x(2,:)))) + exp(1) +20)));
[xk,fk, gradfk_norm, k, xseq, btseq]=newton_bcktrck_verbose(x0, f3, gradf3, Hessf3, kmax, tolgrad, c1, rho, btmax, true);
