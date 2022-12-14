clear all
close all
clc

load("test_functions2.mat")
JFx=fin_diff_J(gradf2, x0, 1e-7, "fw")
