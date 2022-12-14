clear all
close all
clc

load("test_functions2.mat")
h=sqrt(eps);
[gradfx]=fin_diff_grad(f3, x0, h, "c");
[Hessfx]=fin_diff_Hess(f2, x0, sqrt(h));