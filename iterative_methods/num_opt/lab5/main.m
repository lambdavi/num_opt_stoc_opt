clear all
close all
clc

load("test_functions2.mat")
[gradfx]=find_iff_grad(f1, x0, sqrt(eps), "fw");