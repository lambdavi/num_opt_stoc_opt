clear
close all

%c=[0;0];
%r=1;
%x=[-2;0];
%xhat=sphere_projection(x,c,r);
%disp(xhat);

%mins=[1;1];
%maxs=[4;4];
%x=[5;1];
%disp(box_projection(x,mins,maxs));

load("test_functions2.mat")
c1=1e-4;
rho=0.8;
btmax=50;
gamma=0.1;
tolx=1e-12;
c=[-3;-3];
r=1.5;
mins=[-4.5; -4.5];
maxs=[-1.5; -1.5];
handle_sphere = @(x) sphere_projection(x,c,r);
handle_box = @(x) box_projection(x,mins,maxs);

[xk,fk, gradfk_norm, k, xseq, btseq, deltaxk_norm] = untitled( ...
    x0, f2, gradf2, kmax, tolgrad, c1, rho, btmax, gamma, tolx, handle_sphere);


